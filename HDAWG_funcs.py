from time import sleep
from math import ceil
from os import path
import numpy as np
from zhinst.toolkit import Session, CommandTable, Waveforms
import textwrap
import operator

def findfreq(freq,Ecc,loopmax,factor):
    #the single-tone algorithm uses a repeating waveform to reduce the length of waveforms needed to program the AWG.
    #This repeating waveform discretizes the allowed frequency spectrum of the output.
    #Therefore the number of oscillations per repeating segment, and the duration of this segment need to be chosen carefully in order to achieve an accurate frequency.
    #This function accomplishes this using 2 criteria. There is a loopmax, which quantifies the maximum number of samples used for each segment.
    #And Ecc, which quantifies the desired accuracy of the achieved frequency.
    #If the desired accuracy can be reached with less than the maximum number of segments. The minimum number of segments which satisfy the desired accuracy will be chosen.
    #If this first criterion cannot be achieved. The algorithm will choose the parameters for the frequency closest to the desired frequency achievable within the loopmax.
    Tseg=1/(2.4e9)*2**factor
    segprper=1/(freq*Tseg)
    mxloop=int(((loopmax/segprper)))
    

    #Create a interger numpy-array counting up to the maximum number of loops.
    M=np.linspace(1,mxloop,mxloop,dtype=int)

    #First check which must be satisfied is that the segment contains a multiple of 16 samples. This is because of HDAWG playback criteria.
    #If this check is not met, then the segments would be zero-extended to a multiple of 16 samples. 
    #This check also ensures the repeating segments are longer than 100 samples, as in some pathological cases errors might occur otherwise.
    Lf=((M/(freq*Tseg))%16<0.99)*(M/(freq*Tseg)>100)
    Lc=((M/(freq*Tseg))%16>15.01)*(M/(freq*Tseg)>100)
    if True not in Lc and True not in Lf:
        raise(Exception("increase the loopmax"))

    #The second check validates whether each loop number has an achievable frequency close enough to the desired frequency.
    P=np.divide((M/freq%(Tseg)),M)
    Kf=P<1/(freq*Ecc)
    Kc=np.abs(np.divide(Tseg,M)-P)<1/(freq*Ecc)

    #each of these checks produce a truth array, (False,False,True...etc.). Which is multiplied to form two overall truth tables,
    #depending on whether the desired frequency is under the number of samples calculated(ceil) or over the number of samples calculated(floor)
    Kf=Lf*Kf
    Kc=Lc*Kc

    #A&B verify if the first criterion can be reached. 
    A,B= True in Kf, True in Kc
    
    if A and B: #if both ceil and floor can reach first criterion, which can reach it faster?
        Kf=np.where(Kf)
        WN1f=Kf[0][0]+1
        Kc=np.where(Kc)
        WN1c=Kc[0][0]+1
        if WN1c<WN1f:
            C='c'
            WN1=WN1c
        else:
            C='f'
            WN1=WN1f

    #else, just take what works
    elif A:
        K=np.where(Kf)
        WN1=K[0][0]+1
        C='f'

    elif B:
        K=np.where(Kc)
        WN1=K[0][0]+1
        C='c'
    else:#else, use the second criterion.
        Pf=P
        
        Pc=np.abs(np.divide(Tseg,M)-P)
        A=np.min(Pf+np.logical_not(Lf)*10000)
        B=np.min(Pc+np.logical_not(Lc)*10000)
        if A>B:
            WN1=np.where(Pf==np.min(Pf+np.logical_not(Lf)*1000))[0][0]+1
            C='c'
        else:
            WN1=np.where(Pc==np.min(Pc+np.logical_not(Lc)*1000))[0][0]+1
            C='f'
    if C=='f':
        LEN1=int(WN1*2.4e9/(freq*2**factor))
    else:
        LEN1=ceil(WN1*2.4e9/(freq*2**factor))
    
    freqout=WN1/LEN1*24e8/(2**factor)
    return WN1,LEN1,freqout




class Sequence():
    #this class defines and builds the sequences which are loaded onto the HDAWG.
    def __init__(self,Trigchan:int=1,While:bool=True, Ecc:float=1e7, loopmax:int=1e7, acc:float=1e2, t0=0):
        #This initializes the sequence, and defines a few global parameters used for defining the accuracy of the tones, and defines whether the sequence should loop. ("While")
        self.seq = textwrap.dedent(
            """
            while(1){
            waitDigTrigger(DTR);
            }
            """
        )
        
        if While!=True:
            self.seq = self.seq.replace("while(1){","").replace("}","")
            self.While=False
        else:
            self.While=True
        #defines which channel should trigger the sequence.
        self.seq = self.seq.replace("DTR", str(Trigchan)) 
        #defines various values related to the desired frequency accuracy of the HDAWG output.
        self.Ecc = Ecc
        self.lm  = loopmax
        self.acc = acc
        self.t=t0
        #prepares a dictionary which will store the approximated frequencies for the looping segments computed by freqout(). This serves to reduce redundancy.
        self.freqdict={}


    def modulator_on(self):
        #Test code, should turn on modulation on channel 1 and 2.
        if self.While==True:
            self.seq= self.seq[:len(self.seq)-2]

        self.seq+=textwrap.dedent(
            """
            setInt("modulation/0/modulation_mode(0)", 1);
            setInt("modulation/0/modulation_mode(1)", 2);
            
            """
        )
        if self.While==True:
            self.seq +="}"
        
    
    def modulator_off(self):
        #Test code, should turn off modulation on channel 1 and 2.
        if self.While==True:
            self.seq= self.seq[:len(self.seq)-2]

        self.seq+=textwrap.dedent(
            """
            setInt("modulation/0/modulation_mode(0)", 0);
            setInt("modulation/0/modulation_mode(1)", 0);
            
            """
        )
        if self.While==True:
            self.seq +="}"
    
    def addsingletone( #This function produces a single-tone two-channel signal.
        self,
        dur:float       = 0.001,
        freq:float      = 1e6,
        amp:float       = 1.0,
        phase1:float    = None,
        phase2:float    = None,
        Mark:list       = [0,0,0,0],#if and where a marker is desired
        phasetime:float = None, #Time evolution of the phase, this time*freq_real is added onto phase1 and phase2.
        CHAN1:int       = 1, #1st Output channel
        CHAN2:int       = 2  #2nd Output channel
    ):
        #automatically assigns phasetime if it is left undefined
        if phasetime==None:
            phasetime=self.t
        self.t+=dur

        #default phases of channel 1 and 2 are sine and cosine.
        if phase1==None:
            phase1=0
        if phase2==None:
            phase2=phase1+np.pi/2


        #For various calcs, need the following constants.
        per=freq*dur
        factor=max(int(np.log2(2.4e9/(self.acc*freq))),0)
        N=int(round(2.4e9*dur*1/(2**factor),5))

        #If this frequency has already been calculated for, use the previous results, else store the calculation result in the freqdict.
        if freq in self.freqdict:
            [WN1,LEN1,freqout]=self.freqdict[freq]
        else:
            WN1,LEN1,freqout=findfreq(freq,self.Ecc,self.lm,factor)
            self.freqdict[freq]=[WN1,LEN1,freqout]

        #Interpreting part of the output of freqout()
        if WN1==0:
            per=0
        else:
            per=int(per/WN1)
        
        #define phase in terms of a time t0. (This phasetime needs to be assigned. Sequence.sequence_generator takes care of this automatically, but Sequence.addsingletone does not)
        phase1+=(freqout*2*np.pi*phasetime)%(2*np.pi)
        phase2+=(freqout*2*np.pi*phasetime)%(2*np.pi)

        #Find the length of the final, unrepeated, waveform. Checking to ensure the length of this waveform is positive.
        LEN2=N-LEN1*per
        while(LEN2<0):
            per=per-1
            LEN2=N-LEN1*per
        per2=freqout*LEN2/24e8*2**factor

        #Starting the awg_program for this section depending on whether there should be a repeating section:
        if per==0:
            awg_program = textwrap.dedent(
                """
                    playWave(1,marker(LEN2,1)+WAFE1,2,marker(LEN2,1)+WAFE2,3,marker(LEN2,1)+WAFE3,4,marker(LEN2,1)+WAFE4,FACT);
                """
            )
        else:
            awg_program = textwrap.dedent(
                """
                    repeat(PER){
                        playWave(1,marker(LEN1,1)+WAVE1,2,marker(LEN1,1)+WAVE2,3,marker(LEN1,1)+WAVE3,4,marker(LEN1,1)+WAVE4,FACT);
                    }
                    playWave(1,marker(LEN2,1)+WAFE1,2,marker(LEN2,1)+WAFE2,3,marker(LEN2,1)+WAFE3,4,marker(LEN2,1)+WAFE4,FACT);
                """
            )

        
        

        #in order to have arbitrary-channel outputs, the following code makes sure that the outputs: OUT1 & OUT2 end up in the correct spots

        OUT1="+AMP*sine(LEN1, P1, WN1)"
        OUT2="+AMP*sine(LEN1, P2, WN1)"

        
        
        if CHAN1==1:
            awg_program=awg_program.replace("WAVE1","WAVE1"+OUT1)
        elif CHAN1==2:
            awg_program=awg_program.replace("WAVE2","WAVE2"+OUT1)
        elif CHAN1==3:
            awg_program=awg_program.replace("WAVE3","WAVE3"+OUT1)
        else:
            awg_program=awg_program.replace("WAVE4","WAVE4"+OUT1)
        if CHAN2==1:
            awg_program=awg_program.replace("WAVE1","WAVE1"+OUT2)
        if CHAN2==2:
            awg_program=awg_program.replace("WAVE2","WAVE2"+OUT2)
        if CHAN2==3:
            awg_program=awg_program.replace("WAVE3","WAVE3"+OUT2)
        if CHAN2==4:
            awg_program=awg_program.replace("WAVE4","WAVE4"+OUT2)
        awg_program=awg_program.replace("+WAVE1","")
        awg_program=awg_program.replace("+WAVE2","")
        awg_program=awg_program.replace("+WAVE3","")
        awg_program=awg_program.replace("+WAVE4","")
        

        awg_program=awg_program.replace("WAFE","WAVE")
        OUT1="+AMP*sine(LEN2, P1,WN2)"
        OUT2="+AMP*sine(LEN2, P2,WN2)"

        if CHAN1==1:
            awg_program=awg_program.replace("WAVE1","WAVE1"+OUT1)
        elif CHAN1==2:
            awg_program=awg_program.replace("WAVE2","WAVE2"+OUT1)
        elif CHAN1==3:
            awg_program=awg_program.replace("WAVE3","WAVE3"+OUT1)
        else:
            awg_program=awg_program.replace("WAVE4","WAVE4"+OUT1)
        if CHAN2==1:
            awg_program=awg_program.replace("WAVE1","WAVE1"+OUT2)
        if CHAN2==2:
            awg_program=awg_program.replace("WAVE2","WAVE2"+OUT2)
        if CHAN2==3:
            awg_program=awg_program.replace("WAVE3","WAVE3"+OUT2)
        if CHAN2==4:
            awg_program=awg_program.replace("WAVE4","WAVE4"+OUT2)
        awg_program=awg_program.replace("+WAVE1","")
        awg_program=awg_program.replace("+WAVE2","")
        awg_program=awg_program.replace("+WAVE3","")
        awg_program=awg_program.replace("+WAVE4","")
        
        



        #The following code makes sure that the correct markers are retained, while removing all marker waveforms which are unneeded.        
        for i in range(4):
            if Mark[i]==0:
                awg_program=awg_program.replace(str(i+1)+",marker(LEN1,1),","")
                awg_program=awg_program.replace(str(i+1)+",marker(LEN1,1)+",str(i+1)+",")
                awg_program=awg_program.replace(str(i+1)+",marker(LEN2,1),","")
                awg_program=awg_program.replace(str(i+1)+",marker(LEN2,1)+",str(i+1)+",")


        #inserts the correct values for the variable placeholders
        awg_program = awg_program.replace("FACT", str(factor))
        awg_program = awg_program.replace("P1", str(phase1))
        awg_program = awg_program.replace("P2", str(phase2))
        awg_program = awg_program.replace("PER", str(per))
        awg_program = awg_program.replace("WN1", str(WN1))
        awg_program = awg_program.replace("WN2", str(per2))
        awg_program = awg_program.replace("LEN1", str(round(LEN1,5))) #additional rounding is included to mitigate python's floating point weirdness.
        awg_program = awg_program.replace("LEN2", str(round(LEN2,5))) #additional rounding is included to mitigate python's floating point weirdness.
        awg_program = awg_program.replace("AMP", str(amp))
        
        #if the singletone is nested inside of a while loop
        if self.While==True: 
            self.seq= self.seq[:len(self.seq)-2]
        #adding the code, now written, into the sequence
        self.seq+=awg_program

        if self.While==True:
            self.seq +="}"

        
    
    def singletoneIR(
        #This function writes a single-tone channel 1&2 signal+ a channel3 digital modulation 
        # based single-tone into the sequencer. 
        self,
        dur:float       = 0.001,
        freq:float      = 1e6,
        amp:float       = 1.0,
        phase1:float    = None,
        phase2:float    = None,
        IRfreq:float    = None, #frequency of channel 3 output
        IRamp:float     = 1, #amplitude of channel 3 output
        IRphase:float   = 0, #absolute phase of channel 3 output (not mediated by phasetime)
        Mark:list       = [0,0,0,0],#if and where a marker is desired
        phasetime:float = None, #Time evolution of the phase, this time*freq_real is added onto phase1 and phase2.
    ):
        #automatically assigns phasetime if it is left undefined
        if phasetime==None:
            phasetime=self.t
        self.t+=dur
        
        awg_program = textwrap.dedent(#initializes the string written to the sequencer
            """
            """
        )
        if IRfreq!=None: 
            #Checks if the frequency of channel 3 needs to be changed. If it does, a 20 ns break is inserted into the sequence
            # during which the program switches the frequency of the second oscillator, the phase of the 3rd channel,
            # and resets the oscillator phase to 0.
            dur-=2e-7
            awg_program+=textwrap.dedent("""
                playZero(240);
                setDouble("oscs/0/freq", IRfreq);
                setInt("sines/2/phaseshift",P1);
                playZero(240);
                resetOscPhase();""".replace("IRfreq",str(IRfreq)).replace("P1",str(IRphase*180/np.pi))
            )
            phasetime+=2e-7

        

        #default phases of channel 1 and 2 are sine and cosine.        
        if phase1==None:
            phase1=0
        if phase2==None:
            phase2=phase1+np.pi/2


        #For various calcs, need the following factors.
        per=freq*dur
        factor=max(int(np.log2(2.4e9/(self.acc*freq))),0)
        N=int(round(2.4e9*dur*1/(2**factor),5))

        #If this frequency has already been calculated for, use the previous results, else store the calculation result in the freqdict.
        if freq in self.freqdict:
            [WN1,LEN1,freqout]=self.freqdict[freq]
        else:
            WN1,LEN1,freqout=findfreq(freq,self.Ecc,self.lm,factor)
            self.freqdict[freq]=[WN1,LEN1,freqout]

        #Interpreting part of the output of freqout()
        if WN1==0:
            per=0
        else:
            per=int(per/WN1)
        
        #define phase in terms of a time t0. (This phasetime needs to be assigned. Sequence.sequence_generator takes care of this automatically, but Sequence.addsingletone does not)
        phase1+=(freqout*2*np.pi*phasetime)%(2*np.pi)
        phase2+=(freqout*2*np.pi*phasetime)%(2*np.pi)


        #Find the length of the final, unrepeated, waveform. Checking to ensure the length of this waveform is positive.
        LEN2=N-LEN1*per
        while(LEN2<0):
            per=per-1
            LEN2=N-LEN1*per
        per2=freqout*LEN2/24e8*2**factor

        #Starting the awg_program for this section depending on whether there should be a repeating section:
        if per==0:
            awg_program += textwrap.dedent(
                """
                    playWave(1,marker(LEN2,1)+OUTP1,2,marker(LEN2,1)+OUTP2,3,marker(LEN2,1)+IROUTP,4,marker(LEN2,1),FACT);
                """
            )
        else:
            awg_program += textwrap.dedent(
                """
                    repeat(PER){
                        playWave(1,marker(LEN1,1)+OUT1,2,marker(LEN1,1)+OUT2,3,marker(LEN1,1)+IROUTI,4,marker(LEN1,1),FACT);
                    }
                    playWave(1,marker(LEN2,1)+OUTP1,2,marker(LEN2,1)+OUTP2,3,marker(LEN2,1)+IROUTP,4,marker(LEN2,1),FACT);
                """
            )
        
        
        #Ensures the markers are inserted into the correct channels, and placeholders are erased accordingly
        for i in range(4):
            if Mark[i]==0:
                awg_program=awg_program.replace(str(i+1)+",marker(LEN1,1),","")
                awg_program=awg_program.replace(str(i+1)+",marker(LEN1,1)+",str(i+1)+",")
                awg_program=awg_program.replace(str(i+1)+",marker(LEN2,1),","")
                awg_program=awg_program.replace(str(i+1)+",marker(LEN2,1)+",str(i+1)+",")
        
        #The following code replaces the placeholders with the correct waveforms
        OUT1="AMP*sine(LEN1, P1, WN1)"
        OUT2="AMP*sine(LEN1, P2, WN1)"
        IROUTI="rect(LEN1,IRAMP)"
        awg_program=awg_program.replace("OUT1",OUT1)
        awg_program=awg_program.replace("OUT2",OUT2)
        awg_program=awg_program.replace("IROUTI",IROUTI)
        OUTP1="AMP*sine(LEN2, P1,WN2)"
        OUTP2="AMP*sine(LEN2, P2,WN2)"
        IROUTP="rect(LEN2,IRAMP)"
        awg_program=awg_program.replace("OUTP1",OUTP1)
        awg_program=awg_program.replace("OUTP2",OUTP2)
        awg_program=awg_program.replace("IROUTP",IROUTP)

        #inserts the correct values for the variable placeholders
        awg_program=awg_program.replace("IRAMP",str(IRamp))
        awg_program = awg_program.replace("FACT", str(factor))
        awg_program = awg_program.replace("P1", str(phase1))
        awg_program = awg_program.replace("P2", str(phase2))
        awg_program = awg_program.replace("PER", str(per))
        awg_program = awg_program.replace("WN1", str(WN1))
        awg_program = awg_program.replace("WN2", str(per2))
        awg_program = awg_program.replace("LEN1", str(round(LEN1,5)))
        awg_program = awg_program.replace("LEN2", str(round(LEN2,5)))
        awg_program = awg_program.replace("AMP", str(amp))

        
        if self.While==True:#if the singletone is nested inside of a while loop
            self.seq= self.seq[:len(self.seq)-2]
        #adding the code, now written, into the sequence
        self.seq+=awg_program
        if self.While==True:
            self.seq +="}"
        
    def modulatorsingletone(self,tdur=1e-3,freq=None,phase=None,amp=1,Mark=[0,0,0,0]):
        #Writes code into the sequencer for a digital modulator based two-channel single-tone.
        #The code assumes that the digital modulators have been enabled.
        #No method for turning these on/off within the sequencer has been found yet.


        #initial values helpful with further calculations
        per=int(tdur*24e8/240)-1
        LEN1=240
        LEN2=int(tdur*24e8-per*240)
        
        #initializes the program
        awg_program = textwrap.dedent("""
        """)

        #If the frequency is changed, change the frequency.
        if freq!=None:
            if per>0:
                per-=1
            else:
                LEN2-=240
            awg_program+=textwrap.dedent("""
            playZero(176);
            setDouble("oscs/0/freq", %s);
            playZero(64);
            resetOscPhase();""" %str(freq)

            )
        #If the phase is changed, change it accordingly.
        if phase!=None:
            awg_program+=textwrap.dedent("""
            setInt("sines/2/phaseshift",P1);
            setInt("sines/3/phaseshift",P1+90);
            resetOscPhase();""".replace("P1",str(phase*180/np.pi))
            )

        if per>0:#adds the code to produce the rectangular pulse needed for digital modulation, with a check to reduce redundant code
            awg_program += textwrap.dedent(
                """
                    repeat(PER){
                        playWave(1,marker(LEN1,1)+AMP*ones(LEN1),2,marker(LEN1,1)+AMP*ones(LEN1),3,marker(LEN1,1),4,marker(LEN1,1),0);
                    }
                    playWave(1,marker(LEN2,1)+AMP*ones(LEN2),2,marker(LEN2,1)+AMP*ones(LEN2),3,marker(LEN2,1),4,marker(LEN2,1));
                """
            )
        else:
            awg_program += textwrap.dedent(
                """
                    playWave(1,marker(LEN2,1)+AMP*ones(LEN2),2,marker(LEN2,1)+AMP*ones(LEN2),3,marker(LEN2,1),4,marker(LEN2,1));

                """
            )

        #Adds the markers, and removes placeholders
        for i in range(4):
            if Mark[i]==0:
                awg_program=awg_program.replace(str(i+1)+",marker(LEN1,1),","")
                awg_program=awg_program.replace(str(i+1)+",marker(LEN1,1)+",str(i+1)+",")
                awg_program=awg_program.replace(str(i+1)+",marker(LEN2,1),","")
                awg_program=awg_program.replace(str(i+1)+",marker(LEN2,1)+",str(i+1)+",")
                
        #inserts the correct values for the variable placeholders
        awg_program = awg_program.replace("PER", str(per))
        awg_program = awg_program.replace("LEN1", str(round(LEN1,5)))
        awg_program = awg_program.replace("LEN2", str(round(LEN2,5)))
        awg_program = awg_program.replace("AMP", str(amp))

        #if the singletone is nested inside of a while loop
        if self.While==True:
            self.seq=self.seq[:len(self.seq)-2]
        #adding the code, now written, into the sequence
        self.seq+=awg_program
        if self.While==True:
            self.seq+="\n}"

    def addtwotone(
        #Produces a two-tone two-channel signal using the awg-sequencer.
        self,#sequence to append
        dur:float = 0.0001, #duration
        freq1:float = 1e7, #frequency 1st tone
        amp1:float = 0.5, #amplitude 1st tone
        freq2:float = 3e6, #frequency 2nd tone
        amp2:float = 0.5, #amplitude 2nd tone
        phase11:float=0, #phase 1st tone 1st channel
        phase12:float=np.pi/2, #phase 1st tone 2st channel
        phase21:float=0, #phase 2st tone 1st channel
        phase22:float=np.pi/2, #phase 2st tone 2st channel
        Mark:list=[0,0,0,0], #if a marker is desired
        phasetime:float = None, #Time evolution of the phase, this time*freq_real is added onto the phases.
        CHAN1:int       = 1, #1st Output channel
        CHAN2:int       = 2  #2nd Output channel
    ):
        #automatically assigns phasetime if it is left undefined
        if phasetime==None:
            phasetime=self.t
        self.t+=dur

        #factor by which the sampling rate must be reduced (sampling rate=2.4e9/(2**factor))
        factor=max(int(np.log2(2.4e9/(self.acc*max(freq1,freq2)))),0)

        #If this frequency has already been calculated for, use the previous results, else store the calculation result in the freqdict.
        #In the case of the two-tone no approximations are actually needed, but for the sake of consistency, the approximation is still used.
        if freq1 in self.freqdict:
            [_,_,freq1]=self.freqdict[freq1]
        else:
            WN1,LEN1,freq1=findfreq(freq1,self.Ecc,self.lm,factor)
            self.freqdict[freq1]=[WN1,LEN1,freq1]
        if freq2 in self.freqdict:
            [_,_,freq2]=self.freqdict[freq2]
        else:
            WN1,LEN1,freq2=findfreq(freq2,self.Ecc,self.lm,factor)
            self.freqdict[freq2]=[WN1,LEN1,freq2]

    
        #For various calcs, need the following constants.
        per1=freq1*dur
        per2=freq2*dur
        N=N=int(round(2.4e9*dur*1/(2**factor),5))


        #define phase in terms of a time t0. (This phasetime needs to be assigned. Sequence.sequence_generator takes care of this automatically, but Sequence.addtwotone does not)
        phase11+=(freq1*phasetime*np.pi*2)%(2*np.pi)
        phase12+=(freq1*phasetime*np.pi*2)%(2*np.pi)
        phase21+=(freq2*phasetime*np.pi*2)%(2*np.pi)
        phase22+=(freq2*phasetime*np.pi*2)%(2*np.pi)
        
        #initializing code inserted into the sequencer, with all placeholders in place
        SEQ=textwrap.dedent("""
        playWave(1,marker(Num,1)+WAVE1,2,marker(Num,1)+WAVE2,3,marker(Num,1)+WAVE3,4,marker(Num,1)+WAVE4,Fact);""")

        #Replacing the correct placeholders to output in the correct channels.
        OUT1="+add(Amp1*sine(Num,Phase11,Per1),Amp2*sine(Num,Phase21,Per2))"
        OUT2="+add(Amp1*sine(Num,Phase12,Per1),Amp2*sine(Num,Phase22,Per2))"
        if CHAN1==1:
            SEQ=SEQ.replace("WAVE1","WAVE1"+OUT1)
        elif CHAN1==2:
            SEQ=SEQ.replace("WAVE2","WAVE2"+OUT1)
        elif CHAN1==3:
            SEQ=SEQ.replace("WAVE3","WAVE3"+OUT1)
        else:
            SEQ=SEQ.replace("WAVE4","WAVE4"+OUT1)
        if CHAN2==1:
            SEQ=SEQ.replace("WAVE1","WAVE1"+OUT2)
        if CHAN2==2:
            SEQ=SEQ.replace("WAVE2","WAVE2"+OUT2)
        if CHAN2==3:
            SEQ=SEQ.replace("WAVE3","WAVE3"+OUT2)
        if CHAN2==4:
            SEQ=SEQ.replace("WAVE4","WAVE4"+OUT2)
        SEQ=SEQ.replace("+WAVE1","")
        SEQ=SEQ.replace("+WAVE2","")
        SEQ=SEQ.replace("+WAVE3","")
        SEQ=SEQ.replace("+WAVE4","")

        #Ensures the markers are inserted into the correct channels, and placeholders are erased accordingly
        for i in range(4):
            if Mark[i]==0:
                awg_program=awg_program.replace(str(i+1)+",marker(LEN1,1),","")
                awg_program=awg_program.replace(str(i+1)+",marker(Num,1)+",str(i+1)+",")

        #replaces the placeholders with the correct values.
        SEQ=SEQ.replace("Num", str(round(N,5)))
        SEQ=SEQ.replace("Amp1", str(amp1))
        SEQ=SEQ.replace("Amp2", str(amp2))
        SEQ=SEQ.replace("Per1", str(per1))
        SEQ=SEQ.replace("Per2", str(per2))
        SEQ=SEQ.replace("Fact", str(factor))
        SEQ=SEQ.replace("Phase11", str(phase11))
        SEQ=SEQ.replace("Phase12", str(phase12))
        SEQ=SEQ.replace("Phase21", str(phase21))
        SEQ=SEQ.replace("Phase22", str(phase22))


        #if the singletone is nested inside of a while loop
        if self.While==True:
            self.seq=self.seq[:len(self.seq)-2]
        #adding the code, now written, into the sequence
        self.seq+=SEQ
        if self.While==True:
            self.seq+="\n}"

    def addwait(
        #adds a wait to the sequence, with or without markers
        self,
        dur:float = 0.001,#duration of wait
        Mark:list=[0,0,0,0],#if and where a marker is desired
        
    ):
        #adds onto time for phasetime calculations.
        self.t+=dur
        #if the wait is nested inside of a while loop
        if self.While==True:
            self.seq=self.seq[:len(self.seq)-2]
        

        if Mark==[0,0,0,0]: #if no markers are wanted then, the playzero command is used as it takes up close to 0 memory
            self.seq+=textwrap.dedent("""
                repeat(B){
                    playZero(240000000);
                }
            playZero(D);""").replace("D",str(round((dur%0.1)*24e8))).replace("B",str(int(dur*10)))

            
        else:# if there are markers, then there must be a playWave function, which can then be made up of repeated short fragments.
            self.seq+=textwrap.dedent("""
            repeat(NREP){
                playWave(1,marker(240,AUX1),2,marker(240,AUX2),3,marker(240,AUX3),4,marker(240,AUX4));
            }
            playWave(1,marker(N2nd,AUX1),2,marker(N2nd,AUX2),3,marker(N2nd,AUX3),4,marker(N2nd,AUX4));""")

            NREP=int(round(dur*1e7,8))-1
            if NREP==-1:
                NREP=0
                N2nd=int(round(dur*24e8,5))
            else:
                N2nd=int(round(dur*24e8,5))%240+240

            self.seq=self.seq.replace("NREP",str(NREP))
            self.seq=self.seq.replace("N2nd",str(N2nd))
            #A simpler, more memory-intensive method of assigning markers is used here. This is unlikely to cause issues due to the shortness of these fragments.
            self.seq=self.seq.replace("AUX1",str(Mark[0]))
            self.seq=self.seq.replace("AUX2",str(Mark[1]))
            self.seq=self.seq.replace("AUX3",str(Mark[2]))
            self.seq=self.seq.replace("AUX4",str(Mark[3]))
        if self.While==True:
            self.seq+="\n}"
        
        



    
            


    def Upload(self,device,daq):
        #uploads the sequence to the HDAWG.
        #This code is ̶s̶t̶o̶l̶e̶n̶  borrowed from the Zurich Instruments example code. I am not sure how it works, but it has not failed me yet.


        # Create an instance of the AWG Module
        awgModule = daq.awgModule()
        awgModule.set("device", device)
        awgModule.execute()

        # Get the modules data directory
        data_dir = awgModule.getString("directory")
        # All CSV files within the waves directory are automatically recognized by the AWG module
        wave_dir = path.join(data_dir, "awg", "waves")
        if not path.isdir(wave_dir):
            # The data directory is created by the AWG module and should always exist. If this exception
            # is raised, something might be wrong with the file system.
            raise Exception(
                f"AWG module wave directory {wave_dir} does not exist or is not a directory"
            )

        # Save waveform data to CSV
        # Transfer the AWG sequence program. Compilation starts automatically.
        awgModule.set("compiler/sourcestring", self.seq)
        # Note: when using an AWG program from a source file (and only then), the compiler needs to
        # be started explicitly with awgModule.set('compiler/start', 1)
        while awgModule.getInt("compiler/status") == -1:
            sleep(0.1)

        if awgModule.getInt("compiler/status") == 1:
            # compilation failed, raise an exception
            raise Exception(awgModule.getString("compiler/statusstring"))

        if awgModule.getInt("compiler/status") == 0:
            print(
                "Compilation successful with no warnings, will upload the program to the instrument."
            )
        if awgModule.getInt("compiler/status") == 2:
            print(
                "Compilation successful with warnings, will upload the program to the instrument."
            )
            print("Compiler warning: ", awgModule.getString("compiler/statusstring"))
        
        i = 0
        while (awgModule.getDouble("progress") < 1.0) and (
            awgModule.getInt("elf/status") != 1
        ):
            print(f"{i} progress: {awgModule.getDouble('progress'):.2f}")
            #time.sleep(0.2)
            i += 1
        print(f"{i} progress: {awgModule.getDouble('progress'):.2f}")
        if awgModule.getInt("elf/status") == 0:
            print("Upload to the instrument successful.")
        if awgModule.getInt("elf/status") == 1:
            raise Exception("Upload to the instrument failed.")
    


        # This is the preferred method of using the AWG: Run in single mode continuous waveform playback
        # is best achieved by using an infinite loop (e.g., while (true)) in the sequencer program.
        daq.setInt(f"/{device}/awgs/0/single", 1)
        daq.setInt(f"/{device}/awgs/0/enable", 1)


        print("Sequence Activated, Waiting on trigger")
    


    
    def sequence_generator(self,Pulses,Markers):
        Pulses.sort(key=operator.itemgetter("t0"))


        timesp,timesm,tymesp,tymesm=[],[],[],[]
        for P in Pulses:
            timesp.append(P["t0"])
            if P["t0"]<0:
                raise(Exception("Pulses must occur after the trigger"))
            timesm.append(P["t0"]+P["tdur"])
            if P["tdur"]<=0:
                raise(Exception("Pulses must have positive duration"))
        
        
        for M in Markers:
            tymesp.append(M["t0"])
            if M["t0"]<0:
                raise(Exception("Markers must occur after the trigger"))
            tymesm.append(M["t0"]+M["tdur"])
            if M["tdur"]<=0:
                raise(Exception("Pulses must have positive duration"))

        if timesp[0]!=0:
            self.addwait(self,dur=timesp[0],Mark=False)

        for m in range(len(tymesp)):
            if m<len(tymesp)-1:
                if tymesp[m+1]-tymesm[m]<0:
                    raise(Exception("Marker",str(m+1),"and",str(m+2),"overlap. This is not allowed"))

        for t in range(len(timesp)):
            if t<len(timesp)-1:
                if timesp[t+1]-timesm[t]<0:
                    raise(Exception("Pulses",str(t+1),"and",str(t+2),"overlap. This is not allowed"))
        
        Events=[*set(tymesp+tymesm+timesm+timesp)]
        Events.sort()
        Mark=[0,0,0,0]
        
        if Events[0]>0:
            self.addwait(dur=Events[0],Mark=False)
        P=0
        for m in range(len(Events)-1):
            t=Events[m]
            dur=Events[m+1]-t
            if t in tymesp:
                n=tymesp.index(t)
                Mark=Markers[n]["AUX"]
            elif t in tymesm:
                Mark=[0,0,0,0]

            if t in timesp:
                n=timesp.index(t)
                P=Pulses[n]
            elif t in timesm:
                P=0
            
            

            if P==0:
                self.addwait(dur=dur,Mark=Mark)
            else:
                if P["type"]=="Single":
                    freq=P["freq"]
                    P1=P["phase1"]
                    P2=P["phase2"]
                    self.addsingletone(dur=dur,freq=freq,amp=P["Amp"],phase1=P1,phase2=P2,Mark=Mark,phasetime=t)
                elif P["type"]=="Double":
                    freq1=P["freq1"]
                    freq2=P["freq2"]
                    P11=P["phase11"]
                    P12=P["phase12"]
                    P21=P["phase21"]
                    P22=P["phase22"]
                    self.addtwotone(dur=dur,freq1=freq1,amp1=P["Amp1"],freq2=freq2,amp2=P["Amp2"],phase11=P11,phase12=P12,phase21=P21,phase22=P22,Mark=Mark)
            





def add_pulse(
    t0:float=0,                         #Start time in seconds
    tdur:float=2e-4,                    #Duration in seconds
    Type:str="Single",                  #Single/two-tone
    amp1:float=1.0,                     #Amplitude, or Amplitude first tone
    amp2:float=1.0,                     #Amplitude second tone
    freq1:float=1e7,                    #Frequency, or Frequency first tone
    freq2:float=2e7,                    #Frequency second tone
    phase11:float=0,                    #Phase offset first channel first tone
    phase12:float=np.pi/2,              #Phase offset second channel first tone
    phase21:float=0,                    #Phase offset first channel second tone
    phase22:float=np.pi/2,              #Phase offset second channel second tone

):
    P={}
    P["t0"]=t0
    P["tdur"]=tdur
    P["type"]=Type
    if Type=="Single": 
        P["Amp"]=amp1
        P["freq"]=freq1
        P["phase1"]=phase11
        P["phase2"]=phase12
    elif Type=="Double": 
        P["Amp1"]=amp1
        P["freq1"]=freq1
        P["phase11"]=phase11
        P["phase12"]=phase12
        P["Amp2"]=amp2
        P["freq2"]=freq2
        P["phase21"]=phase21
        P["phase22"]=phase22
    else:
        raise(Exception("Type of pulse not Recognized"))
    return P

def add_marker(
    t0:float=0,                         #Start time in seconds best results if 
    tdur:float=2e-4,                    #Duration in seconds
    AUXset:list=[1,0,0,0],
):
    M={}
    M["t0"]=t0
    M["tdur"]=tdur
    M["AUX"]=AUXset
    return M

    