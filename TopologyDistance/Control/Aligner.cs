using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.CompilerServices;
using System.Text;
using System.Threading.Tasks;

namespace TopologyDistance.Control
{
    public class Aligner
    {
        private int[][] SubstitutionMatrix;
        private Dictionary<char, int> SubstitutionMatrixPos;
        public Aligner()
        {
            string PAM30MS = "A R N D C Q E G H I L K M F P S T W Y V B Z J X U *\n" +
                "A 6 -7 -4 -3 -6 -4 -2 -2 -7 -5 -6 -7 -5 -8 -2 0 -1 -13 -8 -2 -7 -6 0 0 0 -17\n" +
                "R -7 8 -6 -10 -8 -2 -9 -9 -2 -5 -7 0 -4 -9 -4 -3 -6 -2 -10 -8 5 -1 0 0 0 -17\n" +
                "N -4 -6 8 2 -11 -3 -2 -3 0 -5 -6 -1 -9 -9 -6 0 -2 -8 -4 -8 -4 -2 0 0 0 -17\n" +
                "D -3 -10 2 8 -14 -2 2 -3 -4 -7 -10 -4 -11 -15 -8 -4 -5 -15 -11 -8 -7 -3 0 0 0 -17\n" +
                "C -6 -8 -11 -14 10 -14 -14 -9 -7 -6 -11 -14 -13 -13 -8 -3 -8 -15 -4 -6 -11 -14 0 0 0 -17\n" +
                "Q -4 -2 -3 -2 -14 8 1 -7 1 -8 -7 -3 -4 -13 -3 -5 -5 -13 -12 -7 -3 4 0 0 0 -17\n" +
                "E -2 -9 -2 2 -14 1 8 -4 -5 -5 -7 -4 -7 -14 -5 -4 -6 -17 -8 -6 -7 -2 0 0 0 -17\n" +
                "G -2 -9 -3 -3 -9 -7 -4 6 -9 -11 -11 -7 -8 -9 -6 -2 -6 -15 -14 -5 -8 -7 0 0 0 -17\n" +
                "H -7 -2 0 -4 -7 1 -5 -9 9 -9 -8 -6 -10 -6 -4 -6 -7 -7 -3 -6 -4 -3 0 0 0 -17\n" +
                "I -5 -5 -5 -7 -6 -8 -5 -11 -9 8 5 -6 -1 -2 -8 -7 -2 -14 -6 2 -6 -7 0 0 0 -17\n" +
                "L -6 -7 -6 -10 -11 -7 -7 -11 -8 5 5 -7 0 -3 -8 -8 -5 -10 -7 0 -7 -7 0 0 0 -17\n" +
                "K -7 0 -1 -4 -14 -3 -4 -7 -6 -6 -7 7 -2 -14 -6 -4 -3 -12 -9 -9 5 4 0 0 0 -17\n" +
                "M -5 -4 -9 -11 -13 -4 -7 -8 -10 -1 0 -2 11 -4 -8 -5 -4 -13 -11 -1 -3 -3 0 0 0 -17\n" +
                "F -8 -9 -9 -15 -13 -13 -14 -9 -6 -2 -3 -14 -4 9 -10 -6 -9 -4 2 -8 -12 -14 0 0 0 -17\n" +
                "P -2 -4 -6 -8 -8 -3 -5 -6 -4 -8 -8 -6 -8 -10 8 -2 -4 -14 -13 -6 -5 -5 0 0 0 -17\n" +
                "S 0 -3 0 -4 -3 -5 -4 -2 -6 -7 -8 -4 -5 -6 -2 6 0 -5 -7 -6 -4 -5 0 0 0 -17\n" +
                "T -1 -6 -2 -5 -8 -5 -6 -6 -7 -2 -5 -3 -4 -9 -4 0 7 -13 -6 -3 -5 -4 0 0 0 -17\n" +
                "W -13 -2 -8 -15 -15 -13 -17 -15 -7 -14 -10 -12 -13 -4 -14 -5 -13 13 -5 -15 -7 -13 0 0 0 -17\n" +
                "Y -8 -10 -4 -11 -4 -12 -8 -14 -3 -6 -7 -9 -11 2 -13 -7 -6 -5 10 -7 -10 -11 0 0 0 -17\n" +
                "V -2 -8 -8 -8 -6 -7 -6 -5 -6 2 0 -9 -1 -8 -6 -6 -3 -15 -7 7 -9 -8 0 0 0 -17\n" +
                "B -7 5 -4 -7 -11 -3 -7 -8 -4 -6 -7 5 -3 -12 -5 -4 -5 -7 -10 -9 5 1 0 0 0 -17\n" +
                "Z -6 -1 -2 -3 -14 4 -2 -7 -3 -7 -7 4 -3 -14 -5 -5 -4 -13 -11 -8 1 4 0 0 0 -17\n" +
                "J 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -17\n" +
                "X 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -17\n" +
                "U 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -17\n" +
                "* -17 -17 -17 -17 -17 -17 -17 -17 -17 -17 -17 -17 -17 -17 -17 -17 -17 -17 -17 -17 -17 -17 -17 -17 1";
            (int[][], Dictionary<char, int>) matrices = LoadSubstitutionMatrixFromString(PAM30MS);
            SubstitutionMatrix = (int[][])matrices.Item1;
            SubstitutionMatrixPos = (Dictionary<char, int>)matrices.Item2;
        }

        public Aligner(string pam30ms)
        {
            (int[][], Dictionary<char, int>) matrices = LoadSubstitutionMatrixFromString(pam30ms);
            SubstitutionMatrix = (int[][])matrices.Item1;
            SubstitutionMatrixPos = (Dictionary<char, int>)matrices.Item2;
        }

        public Aligner(int[][] substitutionMatrix, Dictionary<char, int> substitutionMatrixPos)
        {
            SubstitutionMatrix = substitutionMatrix;
            SubstitutionMatrixPos = substitutionMatrixPos;
        }
        public static int[] GetMaxValue(int[] numbers)
        {
            int maxValue = numbers[0];
            int index = 0;
            int bestIndex = index;
            for (index = 1; index < numbers.Length; index++)
            {
                if (numbers[index] > maxValue)
                {
                    maxValue = numbers[index];
                    bestIndex = index;
                }
            }
            return new int[] { bestIndex, maxValue };
        }
        public (int MaxScorePos, string closestPeptideInProteinDB) GetClosestPeptideInASequence(char[] peptide, char[] protein)
        {
            int[] res = Align(peptide, protein);
            int[] maxScores = GetMaxValue(res);
            int MaxScorePos = maxScores[0];
            int correction = 0;
            if (protein.Length - 1 < MaxScorePos + peptide.Length)
            {
                correction = MaxScorePos + peptide.Length - protein.Length;
            }
            string closestPeptideInProteinDB = new string(protein).Substring(MaxScorePos, MaxScorePos + peptide.Length - correction);
            return (MaxScorePos, closestPeptideInProteinDB);
        }
        public int[] Align(char[] peptide, char[] protein)
        {
            if (protein.Length < peptide.Length)
            {
                Console.WriteLine("Not optimized for peptide larger than proteins");
                return new int[] { 0 };
            }
            int[] alignmentScores = new int[protein.Length];
            for (int x = 0; x < protein.Length; x++)
            {
                int sum = 0;
                for (int y = 0; y < peptide.Length; y++)
                {
                    if (x + y >= protein.Length)
                    {
                        break;
                    }
                    char p1 = peptide[y];
                    char p2 = protein[x + y];
                    int pos1 = SubstitutionMatrixPos[p1];
                    int pos2 = SubstitutionMatrixPos[p2];
                    int add = SubstitutionMatrix[pos1][pos2];
                    sum += add;
                }
                alignmentScores[x] = sum;
            }
            return alignmentScores;
        }
        public static (int[][], Dictionary<char, int>) LoadSubstitutionMatrixFromString(string matrixInText)
        {
            matrixInText = matrixInText.Replace("\r", "");
            int[][] substitutionMatrix = null;
            string[] lines = matrixInText.Split('\n');
            int counter = 0;
            Dictionary<char, int> subsMatrixScores = new Dictionary<char, int>();
            foreach (string line in lines)
            {
                counter++;
                List<string> cols = line.Split(' ').ToList();
                if (counter == 1)
                {
                    List<char> index = new List<char>();
                    foreach (string col in cols)
                    {
                        index.Add(col[0]);
                    }
                    for (int i = 0; i < index.Count; i++)
                    {
                        subsMatrixScores[index[i]] = i;
                    }
                    substitutionMatrix = new int[index.Count + 1][];
                    for (int i = 0; i < substitutionMatrix.Length; i++)
                    {
                        substitutionMatrix[i] = new int[index.Count + 1];
                    }
                }
                else
                {
                    cols.RemoveAt(0);
                    int[] row = new int[cols.Count];
                    for (int i = 0; i < cols.Count; i++)
                    {
                        row[i] = int.Parse(cols[i]);
                    }
                    for (int i = 0; i < row.Length; i++)
                    {
                        substitutionMatrix[i][counter - 2] = row[i];
                    }
                }
            }
            return (substitutionMatrix, subsMatrixScores);
        }
    }
}
