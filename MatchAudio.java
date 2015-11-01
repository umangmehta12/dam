//package matchcalculator;
import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import javax.sound.sampled.*;

public class MatchAudio
{
    public static int FFTLength = 1024;
    public static int HalfFFT = FFTLength/2 ;
    public static int SamplesForOneSecond = 44100;
    public static HashMap<String, SoundFile> PreviousFiles = 
                            new HashMap<String, SoundFile>();

    public static double FrameSize;

    //Purpose: To reomove duplicates from the given array of type double 
    //which is the magnitude from FFTOutput
    //Given: An array of type double which represents the magnitudes
    //Returns: An array of type double after removing duplicates
    public static double[] removeDuplicates(double[] answer)
    {
        double[] cleanedArray = new double[answer.length/2];
        int k = 0;
        for(int i = 0; i<answer.length - (FFTLength);i = i+ FFTLength)
        {
            for(int j = i; j<i+ (FFTLength)/2; j++)
            {
                cleanedArray[k] = answer[j];
                k++;
            }

        }
        System.out.println(k);
        return cleanedArray;
    }
    
    int n, m;
    double[] cos;
    double[] sin;
    double[] window;

    public MatchAudio(int n)
    {
        this.n = n;
        this.m = (int) (Math.log(n) / Math.log(2));
        if (n != (1 << m))
        {
            throw new
            RuntimeException("MatchAudio length must be power of 2");
        }
        cos = new double[n / 2];
        sin = new double[n / 2];
        for (int i = 0; i < n / 2; i++) {
            cos[i] = Math.cos(-2 * Math.PI * i / n);
            sin[i] = Math.sin(-2 * Math.PI * i / n);
        }
        makeWindow();

    }

    protected void makeWindow()
    {
        window = new double[n];
        for (int i = 0; i < window.length; i++)
        {
            window[i] = 0.42 - 0.5 *
            		Math.cos(2 * Math.PI * i / (n - 1))
                    + 0.08 * Math.cos(4 * Math.PI * i / (n - 1));
        }
    }

    public double[] getWindow()
    {
        return window;
    }

    public void MatchAudio(double[] x, double[] y)
    {
        int i, j, k, n1, n2, a;
        double c, s, e, t1, t2;
        j = 0;
        n2 = n / 2;
        for (i = 1; i < n - 1; i++)
        {
            n1 = n2;
            while (j >= n1)
            {
                j = j - n1;
                n1 = n1 / 2;
            }
            j = j + n1;
            if (i < j)
            {
                t1 = x[i];
                x[i] = x[j];
                x[j] = t1;
                t1 = y[i];
                y[i] = y[j];
                y[j] = t1;
            }
        }

        n1 = 0;
        n2 = 1;
        for (i = 0; i < m; i++)
        {

            n1 = n2;
            n2 = n2 + n2;
            a = 0;
            for (j = 0; j < n1; j++)
            {
                c = cos[a];
                s = sin[a];
                a += 1 << (m - i - 1);

                for (k = j; k < n; k = k + n2)
                {
                    t1 = c * x[k + n1] - s * y[k + n1];
                    t2 = s * x[k + n1] + c * y[k + n1];
                    x[k + n1] = x[k] - t1;
                    y[k + n1] = y[k] - t2;
                    x[k] = x[k] + t1;
                    y[k] = y[k] + t2;
                }
            }
        }
    }

    public static void main(String[] args) throws
    IOException, UnsupportedAudioFileException, LineUnavailableException
    {
	//System.out.println("in main");
        String s1 = args[0];
        String s2 = args[1];
        File[] files1 = findOutFiles(s1);
        File[] files2 = findOutFiles(s2);
        dam(files1,files2);
    }

    //Purpose: For each combination of files in both directories calls the 
    //processFiles function
    //Given: Two arrays of type File one for each directory/file 
    public static void dam(File[] files1, File[] files2) 
    throws IOException, UnsupportedAudioFileException, 
           LineUnavailableException
    {
	for(int i = 0; i< files1.length; i++)
        {
            for(int j = 0; j<files2.length; j++)
            {
                SoundFile file1;
                SoundFile file2;
		String file1Name ="/tmp/junk/dir1/".concat(files1[i].getName());
                String file2Name ="/tmp/junk/dir2/".concat(files2[j].getName());
                if(PreviousFiles.containsKey(file1Name))
                {
                    file1 = PreviousFiles.get(file1Name);
		}
                else
                {
		    file1 = new SoundFile(file1Name);
                    PreviousFiles.put(file1Name, file1);
                }
                if(PreviousFiles.containsKey(file2Name))
                {
                    file2 = PreviousFiles.get(file2Name);
		}
                else
                {
                    file2 = new SoundFile(file2Name);
                    PreviousFiles.put(file2Name, file2);
                }
                processFiles(file1, file2);
            }
        }
    }
    
    //Purpose: Finds out all the files in the given directory
    //Given: A String which reprresents the path for a directory
    //Returns: An array of type File
    public static File[] findOutFiles(String s1)
    {
	//System.out.println("in find out files with directory "+ s1);
        File f = new File(s1);
        if(f.isDirectory())
        {
            File[] fArray = f.listFiles();
	    for(int i = 0; i<fArray.length; i++)
		{
			//System.out.println(fArray[i].getName());
		}
            return fArray;
        }
        return null;
    }

    //Given: two SoundFiles
    //Performs a match of the files and returns the corresponding
    //output if its a match or not
    public static void processFiles(SoundFile f1, SoundFile f2)
    {
        MatchAudio MatchAudio = new MatchAudio(FFTLength);
        f1.getPartialMatch(f1,f2);
    }

}
