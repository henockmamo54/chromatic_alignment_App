using Align_Chrom_mzML;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Text.RegularExpressions;
using System.Threading.Tasks;
using XmlMzML;

namespace chromatic_alignment_App
{
    public class Program
    {
        public static void Main(string[] args)
        {
            if (args.Length == 0)
            {
                Console.WriteLine("please provide an input text file.\nYour command should look like:\n> chromatic_alignment_App.exe file.txt\n");
                return;
            }
            else
            {

                var mzmlFileList = readFile(args[0]);

                if (mzmlFileList != null && mzmlFileList.Count > 1)
                    alignChromatograms(mzmlFileList);

            }
        }

        public static List<string> readFile(string path)
        {
            //Extracts the information from files.txt 

            List<string> mzmlFileList = new List<string>();
            if (File.Exists(path))  //check if the file exists 
            {
                Console.WriteLine("reading files.txt");
                try
                {
                    //read all the lines
                    string[] lines = System.IO.File.ReadAllLines(path);
                    lines = lines.Where(x => x.Length > 0).ToArray();

                    foreach (string line in lines)
                    {
                        // remove all extra spaces from the text file. 
                        var temp = line.Trim();

                        if (File.Exists(temp))
                            mzmlFileList.Add(temp);
                        else
                        {
                            Console.WriteLine("Error Reading files.txt \n\t file'" + temp + "' not found.");
                            return null;
                        }

                    }

                }
                catch (Exception e)
                {
                    Console.WriteLine("Error Reading files.txt \n" + e.Message);
                    return null;
                }

            }
            else
            {
                Console.WriteLine("Error Reading files.txt \n\t file not found");
                return null;
            }

            return mzmlFileList;

        }


        public static void alignChromatograms(List<string> mzmlFileList)
        {

            Console.WriteLine("Alignment Started!\n");
            List<Align_Two_Chroms> chromAlignmentPair = new List<Align_Two_Chroms>();
            if (mzmlFileList.Count() < 4)
            {
                batchProcess_mzml_files(Enumerable.Range(0, mzmlFileList.Count()).ToList(), mzmlFileList, chromAlignmentPair);
            }
            else
            {
                int filecount = 0;
                while (filecount < mzmlFileList.Count())
                {
                    var fc_start = filecount == 0 ? 0 : filecount - 1;
                    var fc_end = filecount == 0 ? 3 : 4;

                    if (mzmlFileList.Count() - filecount < 3)
                        fc_end = mzmlFileList.Count() - filecount + 1;

                    var templist = Enumerable.Range(fc_start, fc_end).ToList();
                    //foreach (var kkk in templist)
                    //    Console.WriteLine(kkk.ToString());
                    batchProcess_mzml_files(templist, mzmlFileList, chromAlignmentPair);
                    filecount += 3;
                    //Console.WriteLine("\n");
                }

            }
            Console.WriteLine("Finished Alignment!");

        }
        public static void batchProcess_mzml_files(List<int> index_list, List<string> file_list, List<Align_Two_Chroms> chromAlignmentPair)
        {
            int mz_length = (int)((1500 - 300) * 10 + 1.5);

            bool isTheFirstIndexChromPairRepeated = chromAlignmentPair.Count > 0;

            // prepare chrom pairs
            for (int i = 0; i < index_list.Count - 1; i++)
            {

                //first chrom
                Align_Two_Chroms align = new Align_Two_Chroms();
                if (chromAlignmentPair.Count == 0)
                {
                    align.sFile1 = file_list[index_list[i]];
                    align.mzml1 = new MzML(align.sFile1);
                    align.BuildChromatograms(1, align.mzml1);
                    align.iChrom1_length = align.mzml1.FullScanChromatogram.Count;

                    align.Chrom1 = new float[align.iChrom1_length, mz_length];
                    align.fRet1 = new float[align.iChrom1_length];

                    align.populateChrome(align.mzml1, align.fRet1, align.Chrom1);
                    align.mzml1.CloseFiles();
                    align.mzml1 = null;
                    GC.Collect();
                }
                else
                {
                    align.sFile1 = chromAlignmentPair[chromAlignmentPair.Count - 1].sFile2;
                    align.Chrom1 = chromAlignmentPair[chromAlignmentPair.Count - 1].Chrom2;
                    align.iChrom1_length = chromAlignmentPair[chromAlignmentPair.Count - 1].iChrom2_length;
                    align.fRet1 = chromAlignmentPair[chromAlignmentPair.Count - 1].fRet2;

                }

                // second chrom
                align.sFile2 = file_list[index_list[i + 1]];
                align.mzml2 = new MzML(align.sFile2);
                align.BuildChromatograms(1, align.mzml2);
                align.iChrom2_length = align.mzml2.FullScanChromatogram.Count;

                align.Chrom2 = new float[align.iChrom2_length, mz_length];
                align.fRet2 = new float[align.iChrom2_length];

                align.populateChrome(align.mzml2, align.fRet2, align.Chrom2);
                align.mzml2.CloseFiles();
                align.mzml2 = null;
                GC.Collect();

                chromAlignmentPair.Add(align);
            }

            //execute parallel alignment
            if (isTheFirstIndexChromPairRepeated) chromAlignmentPair.RemoveAt(0);
            run_in_parallel(chromAlignmentPair);

            // clean up the pairs except the last pair. the last pair will be passed to the next batch of pair
            var lastPair = chromAlignmentPair[chromAlignmentPair.Count - 1];
            chromAlignmentPair.Clear();
            chromAlignmentPair.Add(lastPair);


        }

        public static void run_in_parallel(List<Align_Two_Chroms> chromPairs)
        {
            Parallel.ForEach(chromPairs, Align =>
            {

                List<float> result = new List<float>();
                Align.Optimized_Alignment(result);

                //write result to text file
                string[] filename1 = Align.sFile1.Replace(".mzML", "").Split('\\');
                string[] filename2 = Align.sFile2.Replace(".mzML", "").Split('\\');

                //File.WriteAllLines("Align_" + filename1[filename1.Count() - 1] + "_" + filename2[filename2.Count() - 1] + ".txt", result.Select(x => x.ToString()));

                string _filename1 = filename1[filename1.Count() - 1];
                string _filename2 = filename2[filename2.Count() - 1];

                string filecontent = _filename1 + "," + _filename2 + "\n";
                for (int i = 0; i < result.Count - 1; i=i+2)
                {
                    filecontent += result[i] + "," + result[i + 1] + "\n";
                }
                //File.WriteAllLines("Alignment_of_" + filename1[filename1.Count() - 1] + "_" + filename2[filename2.Count() - 1] + ".csv", filecontent);
                using (StreamWriter writer = new StreamWriter("Alignment_of_" + filename1[filename1.Count() - 1] + "_" + filename2[filename2.Count() - 1] + ".csv"))
                {
                    writer.WriteLine(filecontent);
                }
            });
        }

    }
}
