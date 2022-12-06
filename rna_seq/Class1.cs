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
        List<float[]> p;
        List<float[]> c;
        List<float[]> pNorm;
        List<float[]> cNorm;
        List<string> genes;
        const int EXP = 1_000_000;
        const int HISTOGRAM_RANGE = 50;
        const int HISTOGRAM_ROOF = 20;
        const double P_VALUE_FILTER_FLOOR = 0.05d;
        const int MEAN_FILTER_FLOOR = 100;
        const bool IGNORE_ZERO_ROWS = false;
        double[] pVals;
        int pCount;
        int cCount;
        public Class1(string filepath)
        {
            string[] lines = File.ReadAllLines(filepath);
            genes = new();
            pCount = 0;
            cCount = 0;
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
                    p = new List<float[]>();
                    c = new List<float[]>();
                    pNorm = new List<float[]>();
                    cNorm = new List<float[]>();
                }
                else // precist data
                {
                    bool allZeros = true;
                    var tempP = new float[pCount];
                    var tempPNorm = new float[pCount];
                    for (int j = 0; j < pCount; j++)
                    {
                        if (int.Parse(line[j + 2]) > 0)
                        {
                            allZeros = false;
                            tempP[j] = int.Parse(line[j + 2]);
                            tempPNorm[j] = int.Parse(line[j + 2]);
                        }
                    }
                    var tempC = new float[cCount];
                    var tempCNorm = new float[cCount];
                    for (int j = 0; j < cCount; j++)
                    {
                        if (int.Parse(line[j + pCount + 2]) > 0)
                        {
                            allZeros = false;
                            tempC[j] = int.Parse(line[j + pCount + 2]);
                            tempCNorm[j] = int.Parse(line[j + pCount + 2]);
                        }
                    }

                    var meanP = tempP.Sum() / (float)tempP.Length;
                    var meanC = tempC.Sum() / (float)tempC.Length;


                    if (meanP > MEAN_FILTER_FLOOR || meanC > MEAN_FILTER_FLOOR)
                    {
                        p.Add(tempP);
                        pNorm.Add(tempPNorm);
                        c.Add(tempC);
                        cNorm.Add(tempCNorm);
                        genes.Add(line[0]);
                    }
                }
            }

            NormalizeData();
        }
        public void PrintHistogram()
        {
            Histogram h = new Histogram(pVals);


            int[] pocet = new int[HISTOGRAM_RANGE  + 1];

            for (int i = 0; i < pVals.Length; i++)
            {
                int val = (int)Math.Round(pVals[i] * HISTOGRAM_RANGE);
                pocet[val]++;
            }

            int max = pocet.Max();
            int min = pocet.Min();

            var filtered = pocet.Select(p => p).ToArray();
            //var filtered = pocet.ToList().Where(p => p != max).ToArray();

            max = filtered.Max();


            for (int i = 0; i < filtered.Length; i++)
            {
                filtered[i] = (int)Math.Round((filtered[i] - min) / (max - (float)min) * HISTOGRAM_ROOF);
            }

            //Array.Sort(filtered);

            string toPrint = "";
            for (int i = HISTOGRAM_ROOF; i >= 0; i--)
            {
                for (int j = 0; j < filtered.Length; j++)
                {
                    if (filtered[j] >= i)
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
            pVals = new double[p.Count];
            for (int i = 0; i < p.Count; i++)
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
            for (int i = 0; i < pCount; i++)
            {
                xSum += pNorm[index][i];
            }
            for (int i = 0; i < cCount; i++)
            {
                ySum += cNorm[index][i];
            }

            float xAvg = xSum / pCount;
            float yAvg = ySum / cCount;
            GetRozptyl(index, out float sX, out float sY, out float s);

            float sXsq = MathF.Sqrt(sX);
            float sYsq = MathF.Sqrt(sY);


            if (s == 0) s = 1f;

            float t = (xAvg - yAvg) / (s * MathF.Sqrt((1f / pCount) + (1f / cCount)));

            return t;
        }
        private double GetPvalue(int index)
        {

            double t = GetTstat(index);

            var pp = new List<double>();
            for (int i = 0; i < pCount; i++)
            {
                pp.Add(pNorm[index][i]);
            }
            var cc = new List<double>();
            for (int i = 0; i < cCount; i++)
            {
                cc.Add(cNorm[index][i]);
            }

            var test = new TTest(t, pCount + cCount - 2);
            if (test.PValue >= 1) { }
            return test.PValue > 1 ? 1 : test.PValue;
        }
        private float GetFoldChange(int index)
        {
            return MathF.Log2(pCount / (float)cCount);
        }
        private void GetRozptyl(int index, out float sX, out float sY, out float s)
        {
            sX = 0;
            sY = 0;

            float xSum = 0;
            float ySum = 0;
            for (int i = 0; i < pCount; i++)
            {
                xSum += pNorm[index][i];
            }
            for (int i = 0; i < cCount; i++)
            {
                ySum += cNorm[index][i];
            }

            float xAvg = xSum / pCount;
            float yAvg = ySum / cCount;


            for (int i = 0; i < pCount; i++)
            {
                sX += (pNorm[index][i] - xAvg) * (pNorm[index][i] - xAvg);
            }
            for (int i = 0; i < cCount; i++)
            {
                sY += (cNorm[index][i] - yAvg) * (cNorm[index][i] - yAvg);
            }

            sX /= pCount;
            sY /= cCount;

            float sPwr = 0;

            sPwr = ((pCount - 1)*sX + (cCount - 1)*sY) / (pCount + cCount - 2);

            s = MathF.Sqrt(sPwr);
        }
        private void NormalizeData()
        {
            for (int i = 0; i < pCount; i++)
            {
                float sum = 0;
                for (int j = 0; j < p.Count; j++)
                {
                    sum += p[j][i];
                }
                for (int j = 0; j < p.Count; j++)
                {
                    pNorm[j][i] = p[j][i] / sum * EXP;
                }
            }
            for (int i = 0; i < cCount; i++)
            {
                float sum = 0;
                for (int j = 0; j < c.Count; j++)
                {
                    sum += c[j][i];
                }
                for (int j = 0; j < c.Count; j++)
                {
                    cNorm[j][i] = c[j][i] / sum * EXP;
                }
            }
        }

        public void Filter()
        {
            var pvalss = pVals.ToList();
            var geness = genes.ToList();


            for (int i = pVals.Length - 1; i >= 0; i--)
            {
                if (pVals[i] < P_VALUE_FILTER_FLOOR)
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
                    if (pvalss[j] > pvalss[j + 1])
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
