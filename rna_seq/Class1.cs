using System;
using System.Collections.Generic;
using System.Diagnostics.Metrics;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Accord.Math;
using Accord.Statistics.Testing;
using Accord.Statistics.Visualizations;

namespace rna_seq
{
    internal class Class1
    {
        float[,] p;
        float[,] c;
        float[,] pNorm;
        float[,] cNorm;
        string[] genes;
        const int EXP = 1_000_000;
        double[] pVals;
        public Class1(string filepath)
        {
            string[] lines = File.ReadAllLines(filepath);
            genes = new string[lines.Length - 1];
            int pCount = 0;
            int cCount = 0;
            for (int i = 0; i < lines.Length; i++)
            {
                var line = lines[i].Split(';');
                
                //spocitat pocet P a C
                if (i == 0)
                {
                    foreach (var column in line)
                    {
                        if (column.StartsWith('P')) pCount++;
                        if (column.StartsWith('C')) cCount++;
                    }
                    p = new float[lines.Length - 1, pCount];
                    c = new float[lines.Length - 1, cCount];
                    pNorm = new float[lines.Length - 1, pCount];
                    cNorm = new float[lines.Length - 1, cCount];
                }
                else // precist data
                {
                    genes[i - 1] = line[0];
                    for (int j = 0; j < pCount; j++)
                    {
                        p[i-1, j] = int.Parse(line[j + 2]);
                        pNorm[i-1, j] = int.Parse(line[j + 2]);
                    }
                    for (int j = 0; j < cCount; j++)
                    {
                        c[i - 1, j] = int.Parse(line[j + pCount + 2]);
                        cNorm[i - 1, j] = int.Parse(line[j+pCount + 2]);
                    }
                }
            }

            NormalizeData();
        }
        public void PrintHistogram()
        {
            Histogram h = new Histogram(pVals);

            const int histogramRange = 20;

            int[] pocet = new int[histogramRange  + 1];

            for (int i = 0; i < pVals.Length; i++)
            {
                int val = (int)Math.Round(pVals[i] * histogramRange);
                pocet[val]++;
            }

            for (int i = 0; i < pVals.Length; i++)
            {
                int val = (int)Math.Round(pVals[i] * histogramRange);
                if (val == 0) { }
                pocet[val]++;
            }

            int max = pocet.Max();
            int min = pocet.Min();

            //var filtered = pocet.Select(p => p).ToArray();
            var filtered = pocet.ToList().Where(p => p != max).ToArray();

            max = filtered.Max();


            for (int i = 0; i < filtered.Length; i++)
            {
                filtered[i] = (int)Math.Round((filtered[i] - min) / (max - (float)min) * histogramRange);
            }

            Array.Sort(filtered);

            string toPrint = "";
            for (int j = filtered.Length - 1; j >= 0; j--)
            {
                for (int i = 0; i < histogramRange; i++)
                {
                    if (filtered[i] >= j)
                    {
                        toPrint += '|';
                    }
                    else
                    {
                        toPrint += ' ';
                    }
                }
                toPrint += '\n';
            }

            Console.WriteLine(toPrint);
        }
        public void Ttest()
        {
            pVals = new double[p.GetLength(0)];
            for (int i = 0; i < p.GetLength(0); i++)
            {
                double pVal = GetPvalue(i);
                double fCh = GetFoldChange(i);

                pVals[i] = pVal;
            }
        }
        private double GetTstat(int index)
        {

            float xSum = 0;
            float ySum = 0;
            for (int i = 0; i < p.GetLength(1); i++)
            {
                xSum += pNorm[index, i];
            }
            for (int i = 0; i < c.GetLength(1); i++)
            {
                ySum += cNorm[index, i];
            }

            float xAvg = xSum / p.GetLength(1);
            float yAvg = ySum / c.GetLength(1);
           
            GetRozptyl(index, out float sX, out float sY, out float s);

            float sXsq = MathF.Sqrt(sX);
            float sYsq = MathF.Sqrt(sY);


            if (s == 0) s = 1f;

            float t = (xAvg - yAvg) / (s * MathF.Sqrt((1f / p.GetLength(1)) + (1f / c.GetLength(1))));

            return t;
        }
        private double GetPvalue(int index)
        {

            double t = GetTstat(index);

            var pp = new List<double>();
            for (int i = 0; i < p.GetLength(1); i++)
            {
                pp.Add(pNorm[index, i]);
            }
            var cc = new List<double>();
            for (int i = 0; i < c.GetLength(1); i++)
            {
                cc.Add(cNorm[index, i]);
            }

            var test = new TTest(t, p.GetLength(1) + c.GetLength(1) - 2);
            if (test.PValue >= 1) { }
            return test.PValue > 1 ? 1 : test.PValue;
        }
        private float GetFoldChange(int index)
        {
            return MathF.Log2(p.GetLength(1) / (float)c.GetLength(1));
        }
        private void GetRozptyl(int index, out float sX, out float sY, out float s)
        {
            sX = 0;
            sY = 0;

            float xSum = 0;
            float ySum = 0;
            for (int i = 0; i < p.GetLength(1); i++)
            {
                xSum += pNorm[index, i];
            }
            for (int i = 0; i < c.GetLength(1); i++)
            {
                ySum += cNorm[index, i];
            }

            float xAvg = xSum / p.GetLength(1);
            float yAvg = ySum / c.GetLength(1);


            for (int i = 0; i < p.GetLength(1); i++)
            {
                sX += (pNorm[index, i] - xAvg) * (pNorm[index, i] - xAvg);
            }
            for (int i = 0; i < c.GetLength(1); i++)
            {
                sY += (cNorm[index, i] - yAvg) * (cNorm[index, i] - yAvg);
            }

            sX /= p.GetLength(1);
            sY /= c.GetLength(1);

            float sPwr = 0;

            sPwr = ((p.GetLength(1) - 1)*sX + (c.GetLength(1) - 1)*sY) / (p.GetLength(1) + c.GetLength(1) - 2);

            s = MathF.Sqrt(sPwr);
        }
        private void NormalizeData()
        {
            pNorm = new float[p.GetLength(0), p.GetLength(1)];
            cNorm = new float[c.GetLength(0), c.GetLength(1)];

            for (int i = 0; i < p.GetLength(1); i++)
            {
                float sum = 0;
                for (int j = 0; j < p.GetLength(0); j++)
                {
                    sum += p[j, i];
                }
                for (int j = 0; j < p.GetLength(0); j++)
                {
                    pNorm[j, i] = p[j, i] / sum * EXP;
                }
            }
            for (int i = 0; i < c.GetLength(1); i++)
            {
                float sum = 0;
                for (int j = 0; j < c.GetLength(0); j++)
                {
                    sum += c[j, i];
                }
                for (int j = 0; j < c.GetLength(0); j++)
                {
                    cNorm[j, i] = c[j, i] / sum * EXP;
                }
            }
        }

        public void Filter()
        {
            var pvalss = pVals.ToList();
            var geness = genes.ToList();


            for (int i = pVals.Length - 1; i >= 0; i--)
            {
                if (pVals[i] < 0.05d)
                {
                    pvalss.RemoveAt(i);
                    geness.RemoveAt(i);
                }
            }
            bool swapRequired;
            for (int i = 0; i < pvalss.Count - 1; i++)
            {
                swapRequired = false;
                for (int j = 0; j < pvalss.Count - i - 1; j++)
                    if (pvalss[j] < pvalss[j + 1])
                    {
                        var tempVar = pvalss[j];
                        pvalss[j] = pvalss[j + 1];
                        pvalss[j + 1] = tempVar;
                        swapRequired = true;
                        var tempVarr = geness[j];
                        geness[j] = geness[j + 1];
                        geness[j + 1] = tempVarr;
                    }
                if (swapRequired == false)
                    break;
            }

            for (int i = 0; i < 10; i++)
            {
                Console.WriteLine(geness[i]);
            }

        }
    }
}
