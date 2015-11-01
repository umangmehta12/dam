//package match;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.util.ArrayList;
import java.util.HashMap;
import javax.sound.sampled.AudioFormat;
import javax.sound.sampled.AudioInputStream;
import javax.sound.sampled.AudioSystem;
import javax.sound.sampled.LineUnavailableException;
import javax.sound.sampled.UnsupportedAudioFileException;

public class SoundFile 
{
    public static int FFTLength = 1024;
    public static int HalfFFT = FFTLength/2 ;
    public static int ChunksForHalfSecond = 22;
    
    public SoundFile(String s) throws IOException
            ,UnsupportedAudioFileException, LineUnavailableException
    {
	try
	{
        this.FilePath = s;
        //this.File = new File(s);
        this.FileArray =  getbytes(s);
        //this.Duration = getDuration();
        //double[] imaginary1 = getInitialImaginaryArray(this.FileArray.length);
        this.MagnitudeArray = getFFTOutput(this.FileArray);
        this.FingerPrint = constructFingerPrint(this.MagnitudeArray);
        this.NameToDisplay = this.getTrimmedName();
        this.FingerHash = getFingerHash(this);
        if(s.contains("_mp3"))
        {
            String tempName = this.NameToDisplay;
            this.NameToDisplay = tempName.replace("_mp3.wav", ".mp3");
        }
        else if(s.contains("_ogg"))
        {
            String tempName = this.NameToDisplay;
            this.NameToDisplay = tempName.replace("_ogg.wav", ".ogg");
        }
	}
	catch(Exception e)
	{
		System.out.println(s);
		e.printStackTrace();
	}
    }
    
    public String FilePath;
    public String NameToDisplay;
    public double Duration;
    public HashMap<Integer,Peak> FingerPrint;
    public double[] FileArray;
    public double[] MagnitudeArray;
    public File File;  
    public HashMap<Integer,Peak>FingerHash;

    //Purpose: To get duration of the input file
    //Given: A file
    //Returns: Double which specifies the duration of the file
    public double getDuration()throws 
           UnsupportedAudioFileException, IOException, LineUnavailableException
    {
        File file = new File(this.FilePath);
        AudioInputStream audioInputStream = 
               AudioSystem.getAudioInputStream(file);
       AudioFormat format = audioInputStream.getFormat();
       float rate = format.getFrameRate();
       long frames = audioInputStream.getFrameLength();
       double durationInSeconds = frames / rate ;
       return durationInSeconds;
   }
   
    //Purpose: To get a hashmap which serves as the fingerprint of the chunk
    //Given: Array of type double which contains the magnitudes
    //Returns: A hashmap<Integer, Peak> which stores ChunkId as the key and 
    //corresponding Peak as the value
    public static HashMap<Integer, Peak>
    constructFingerPrint(double[] magnitudeArray)
    {
        HashMap<Integer,Peak> FingerPrint = new HashMap<Integer, Peak>();
        int k = 1;
        for(int i=0; i<=magnitudeArray.length - HalfFFT; i = i+HalfFFT)
        {
            double peak = 0.0;
            int peakIndex = 0;
            for(int j=i; j< i+HalfFFT;j++)
            {
                if(magnitudeArray[j] >= peak)
                {
                    peak = magnitudeArray[j];
                    peakIndex = j;
                    continue;
                }
            }

            Peak P = new Peak();
            P.Index = (peakIndex - ((k-1)*512));
            P.Value = peak;
            FingerPrint.put(k, P);
            k++;
        }

        return FingerPrint;
    }
    
    //Purpose: To get a magnitude array from the real and imaginary parts
    //after performing FFT and removing unwanted frequencies
    //Given: Two arrays for real and imaginary parts both of type double
    //Returns: An array of type double represents magnitudes after removing 
    //unwanted frequencies.
    public static double[]
    		getFFTOutput(double[] doublearray1)
    {
        MatchAudio MatchAudio = new MatchAudio(FFTLength);
        double[] answer = new double[doublearray1.length/2];

        int ctr = 0;
        for(int i = 0; i<=doublearray1.length - FFTLength; i = i+ FFTLength)
        {
            int k = 0;
            double[] fftReal = new double[1024];
            double[] fftImg = new double[1024];
            double[] magArray = new double[512];
            for(int j = i; j< i+FFTLength; j++)
            {
                fftReal[k] = doublearray1[j];
                fftImg[k] = 0;
                k++;
            }
            beforeAfter(MatchAudio, fftReal, fftImg);
            magArray = getMagnitudeArray(fftReal, fftImg);

            for(int l = 0; l<magArray.length; l++)
            {
                answer[ctr] = magArray[l];
                ctr++;
            }
        }

        return answer;
    }

    //Given: real and imaginary part of array both of type double 
    //Returns: magnitude array of type double
    public static double[] getMagnitudeArray
    (double[] ubytearray1, double[] im1)
    {
        //int j = 0;
        double[] magnitudeArray = new double[ubytearray1.length/2];
        for(int i = 0; i < (HalfFFT - 1); i++)
        {
            if(i <= 25 || i > 465)
            {
                magnitudeArray[i] = 0;
                //j++;
                continue;
            }
            double magnitude = getMagnitude(ubytearray1[i], im1[i]);
            magnitudeArray[i] = magnitude;
            //j++;
        }

       //double[] cleanedArray = removeUnwatedEntries(magnitudeArray);
       //return cleanedArray;
        return magnitudeArray;
    }
    
    //Given: Magnitude array of type double
    //Returns: Magnitude array of type double after removing duplicates 
    //and frequencies those are not audible to human ears
    public static double[] removeUnwatedEntries(double[] magnitudeArray)
    {
        int k = 0;
        double[] answer = new double[magnitudeArray.length/2];

            for(int j = 0; j<= (HalfFFT - 1); j++)
            {
                if(j<=5)
                {
                    answer[k] = 0;
                    k++;
                    continue;
                }
                if(j>465) //20,000 Hz
                {
                    answer[k] = 0;
                    k++;
                    continue;
                }
                answer[k] = magnitudeArray[j];
                k++;
            }
        return answer;
    }
    
    //Given: Real and imaginary parts of a number both of type double
    //Returns: Magnitude of type double calculated from the given real 
    //and imaginary data
    public static double getMagnitude(double ubytearray1, double im1)
    {
        return Math.sqrt(Math.pow(ubytearray1, 2) + Math.pow(im1, 2));
    }
    
    //Given: Length as integer
    //Returns: array of double that represents imaginary array        
    public static double[] getInitialImaginaryArray(int length)
    {
       double[] answer = new double[length];
       for(int i = 0; i< length; i++)
       {
           answer[i] = 0.0;
       }
       return answer;
    }
    
    //Purpose: Calls the MatchAudio constructor to
    //perform match of audio
    //Given: Instance of MatchAudio, two arrays of type double representing
    //real and imaginary parts
    public static void beforeAfter(MatchAudio MatchAudio,
    		double[] re, double[] im)
    {
        MatchAudio.MatchAudio(re, im);
    }
    
    //Purpose: To trim the given input
    //Given: A String    
    //Returns: returns a trimmed string
    public String getTrimmedName() 
    {
        String s = this.FilePath;
        String[] tempArray;
        if (s.contains("/")) {
            tempArray = s.split("/");
            s = getTrimmedNameHelper(tempArray);
            return s;
        }

        return s;
    }
    
    //Purpose: Gets the name of the file from the given input
    //Given: Array of strings
    //Returns: string (which is name of the file)
    public String getTrimmedNameHelper(String[] tempArray) 
    {
        for (int i = 0; i <= tempArray.length; i++) {
            if (tempArray[i].contains(".wav")) {
                return tempArray[i];
            }
        }
        return " ";
    }
    
    //Purpose: To get the byte representation of the file
    //Given: A file
    //Returns: array of type double which has the byte
    //representation for the input file
    public static double[] getbytes(String s) throws IOException
    {
        File file = new File(s);
        FileInputStream fileInputStream = null;
        byte[] bFile = new byte[(int) file.length()];
        
	byte[] bFileDouble = new byte[(int) file.length() - 44];
	
        fileInputStream = new FileInputStream(file);
        fileInputStream.read(bFile);
        fileInputStream.close();

        int j = 0;
        for (int i = 44; i < bFile.length; i++) 
        {
            bFileDouble[j] = ((bFile[i]));
            j++;
        }
        
        double[] doubleArray = getdouble(bFileDouble);
        return doubleArray;
    }
    
    //Given: A byte array
    //Returns: An array of type double which is the representation by 
    //wrapping the array contents using the little endian convention
    public static double[] getdouble(byte[] bArray)
    {
        double[] doubleArray = new double[(int)(bArray.length/2)];
        
        int j = 0;
        for(int i=0;i<bArray.length - 2;i+=2)
        { 
            byte[] arr = new byte[2];
            arr[0] = bArray[i];
            arr[1] = bArray[i+1];
            double val;
            val= 
            ByteBuffer.wrap(arr).order(ByteOrder.LITTLE_ENDIAN).getShort();
            doubleArray[j] = val;
            j++;
            
        }
        return doubleArray;
    }

    //Purpose: Calculates a custom partial match algorithm, also
    //replaces the .mp3 files converted to .wav to their
    //corresponding format for the required output of the
    //program
    //Given: Two SoundFile instances of the files to be compared
//    public void getPartialMatch(SoundFile f1,SoundFile f2)
//    {
//        HashMap<Integer, Peak> FH1 = f1.FingerHash;
//        HashMap<Integer, Peak> FH2 = f2.FingerHash;
//	
//

//
//        for(int i = 0; i<(FH1.size()) - 20; i = i+1)
//        {
//            Peak p1 = null;
//            Peak p2 = null;
//
//            for(int j = 0; j<(FH2.size()) - 20; j = j+1)
//            {
//		//System.out.println("i " + i);
//		//System.out.println("j " + j);
//                p1 = FH1.get(i);
//                p2 = FH2.get(j);
//                if(p1.Index == p2.Index)
//                {
//                    if((FH1.get(i+1).Index == FH2.get(j+1).Index)&&
//			(FH1.get(i+2).Index == FH2.get(j+2).Index)&&
//                        (FH1.get(i+3).Index == FH2.get(j+3).Index)&&
//			(FH1.get(i+4).Index == FH2.get(j+4).Index)&&
//			(FH1.get(i+5).Index == FH2.get(j+5).Index)&&
//			(FH1.get(i+6).Index == FH2.get(j+6).Index)&&
//                        (FH1.get(i+7).Index == FH2.get(j+7).Index)&&
//                        (FH1.get(i+9).Index == FH2.get(j+9).Index)&&
//                        (FH1.get(i+11).Index == FH2.get(j+11).Index)&&
//                        (FH1.get(i+13).Index == FH2.get(j+13).Index)&&
//			(FH1.get(i+15).Index == FH2.get(j+15).Index)&&
//			(FH1.get(i+16).Index == FH2.get(j+16).Index)&&
//                        (FH1.get(i+17).Index == FH2.get(j+17).Index)&&
//			(FH1.get(i+18).Index == FH2.get(j+18).Index)&&
//                        (FH1.get(i+19).Index == FH2.get(j+19).Index) 
//                            )
//                    {
//                        printMatch(f1,f2,i,j);
//                        return;
//                    }
//                }
//            }
//        }
//    }
    
    public void getPartialMatch(SoundFile f1,SoundFile f2)
    {
        HashMap<Integer, Peak> FH1 = f1.FingerHash;
        HashMap<Integer, Peak> FH2 = f2.FingerHash;
	for(int i = 0; i<(FH1.size()) - 20; i = i+1)
        {
            Peak p1 = null;
            Peak p2 = null;

            for(int j = 0; j<(FH2.size()) - 20; j = j+1)
            {
		//System.out.println("i " + i);
		//System.out.println("j " + j);
                p1 = FH1.get(i);
                p2 = FH2.get(j);
                int count = 0;
                if(p1.Index == p2.Index)
                {
                    if(Math.abs(FH1.get(i+1).Index - FH2.get(j+1).Index) <= 0)
                    {
                        count++;
                    }
                    if(Math.abs(FH1.get(i+2).Index - FH2.get(j+2).Index) <= 0)
                    {
                        count++;
                    }
                    if(Math.abs(FH1.get(i+3).Index - FH2.get(j+3).Index) <= 0)
                    {
                        count++;
                    }
                    if(Math.abs(FH1.get(i+4).Index - FH2.get(j+4).Index) <= 0)
                    {
                        count++;
                    }
                    if(Math.abs(FH1.get(i+5).Index - FH2.get(j+5).Index) <= 0)
                    {
                        count++;
                    }
                    if(Math.abs(FH1.get(i+6).Index - FH2.get(j+6).Index) <= 0)
                    {
                        count++;
                    }
                    if(Math.abs(FH1.get(i+7).Index - FH2.get(j+7).Index) <= 0)
                    {
                        count++;
                    }
                    if(Math.abs(FH1.get(i+8).Index - FH2.get(j+8).Index) <= 0)
                    {
                        count++;
                    }
                    if(Math.abs(FH1.get(i+9).Index - FH2.get(j+9).Index) <= 0)
                    {
                        count++;
                    }
                    if(Math.abs(FH1.get(i+10).Index - FH2.get(j+10).Index) <= 0)
                    {
                        count++;
                    }
                    if(Math.abs(FH1.get(i+11).Index - FH2.get(j+11).Index) <= 0)
                    {
                        count++;
                    }
                    if(Math.abs(FH1.get(i+12).Index - FH2.get(j+12).Index) <= 0)
                    {
                        count++;
                    }
                    if(Math.abs(FH1.get(i+13).Index - FH2.get(j+13).Index) <= 0)
                    {
                        count++;
                    }
                    if(Math.abs(FH1.get(i+14).Index - FH2.get(j+14).Index) <= 0)
                    {
                        count++;
                    }
                    if(Math.abs(FH1.get(i+15).Index - FH2.get(j+15).Index) <= 0)
                    {
                        count++;
                    }
                    if(Math.abs(FH1.get(i+16).Index - FH2.get(j+16).Index) <= 0)
                    {
                        count++;
                    }
                    if(Math.abs(FH1.get(i+17).Index - FH2.get(j+17).Index) <= 0)
                    {
                        count++;
                    }
                    if(Math.abs(FH1.get(i+18).Index - FH2.get(j+18).Index) <= 0)
                    {
                        count++;
                    }
                    if(Math.abs(FH1.get(i+19).Index - FH2.get(j+19).Index) <= 0)
                    {
                        count++;
                    }
                    
                    if(count >14)
                    {
                        printMatch(f1,f2,i,j);
                        return;
                    }
                }
            }
        }
    }
    
//    public void getPartialMatch(SoundFile f1,SoundFile f2)
//    {
//        HashMap<Double, Peak> FH1 = getFingerHash(f1);
//        HashMap<Double, Peak> FH2 = getFingerHash(f2);
//
//        for(double i = 0; i<(FH1.size()) - 10; i = i+1)
//        {
//            Peak p1 = null;
//            Peak p2 = null;
//
//            for(double j = 0; j<(FH2.size()) - 10; j = j+1)
//            {
//                p1 = FH1.get(i);
//                p2 = FH2.get(j);
//                if(p1.Index == p2.Index)
//                {
//                    if((FH1.get(i+1).Index == FH2.get(j+1).Index)&&
//                        (FH1.get(i+3).Index == FH2.get(j+3).Index)&&
//			(FH1.get(i+5).Index == FH2.get(j+5).Index)&&
//                        (FH1.get(i+7).Index == FH2.get(j+7).Index)&&
//                        (FH1.get(i+9).Index == FH2.get(j+9).Index))
//                    {
//                        printMatch(f1,f2,i,j);
//                        return;
//                    }
//                }
//            }
//        }
//    }
    
    //Purpose: To print match if two audio fragments sound similar
    //Given: Two SoundFiles and two double which represents time offsets
    //for the fragments which matches in both files
    public static void printMatch(SoundFile f1,SoundFile f2,double i,double j)
    {	
	double t1 = (double)i/4;
	double t2 = (double)j/4;
        System.out.println("MATCH " + f1.NameToDisplay+ " "+ f2.NameToDisplay+
                " "+ t1 +" "+t2);
    }
    
    //Purpose: To get a summarized fingerprint
    //Given: A SoundFile
    //Returns: A HashMap<Double, Peak> where key is the time and value is Peak
    public HashMap<Integer,Peak> getFingerHash(SoundFile f1)
    {
	//System.out.println("hello");
        int time = 0;
        HashMap<Integer, Peak> FingerHash = new HashMap<Integer, Peak>();
        for(int i = 1; i<f1.FingerPrint.size() - (ChunksForHalfSecond/2); 
                i = (i+ChunksForHalfSecond/2))
        {
            double max = 0;
            Peak maxPeak = null;
            for(int j = i; j<i+(ChunksForHalfSecond/2); j++)
            {
		
                Peak fp = f1.FingerPrint.get(j);
                if(fp.Value > max)
                {
                    max = fp.Value;
                    maxPeak = fp;
                }
            }
            if(maxPeak == null)
            {
                Peak dummyPeak = new Peak();
                dummyPeak.Index = 0;
                dummyPeak.ChunkId = 0;
                dummyPeak.Value = 0;
                FingerHash.put(time, dummyPeak);
            }
            else
            {
                FingerHash.put(time, maxPeak);
	    }
		
	
            time = time+1;
        }
        return FingerHash;
    }
}
