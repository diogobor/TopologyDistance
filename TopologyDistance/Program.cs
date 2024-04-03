using System.Runtime.InteropServices;
using System.Text;
using System.Text.RegularExpressions;
using TopologyDistance.Control;
using TopologyDistance.Model;
using static System.Formats.Asn1.AsnWriter;

namespace TopologyDistance
{
    internal class Program
    {
        private static string[] HeaderLineCSV { get; set; }
        static void Main(string[] args)
        {
            List<PPI> ppis = new();
            //string filePath = @"D:\Experimental_data\Cong\cytoscapedefaultnode.csv";
            string filePath = @"D:\Experimental_data\Cong\cytoscape.csv";
            string fastaFile = @"D:\Experimental_data\Cong\20220323_CW_DSBSO_Paper_HEK_db.fasta";

            #region read fasta
            string[] fasta_lines = File.ReadAllLines(fastaFile);

            List<FastaItem> fasta_items = new List<FastaItem>();
            string accession_number = "";
            bool canAdd = false;

            StringBuilder sb_seq = new();
            foreach (var fasta_line in fasta_lines)
            {
                if (fasta_line.StartsWith(">"))
                {
                    if (canAdd)
                    {
                        string sequence = sb_seq.ToString().Replace("\n", "").Replace("\r", "").Trim();
                        FastaItem _fasta = new(accession_number.Trim(), sequence, sequence.Length);
                        sb_seq.Clear();
                        fasta_items.Add(_fasta);
                    }
                    accession_number = Regex.Split(fasta_line, " ")[0].Replace(">", "");
                    if (accession_number.Contains("|"))
                        accession_number = Regex.Split(accession_number, "\\|")[1];

                    canAdd = true;
                }
                else
                    sb_seq.Append(fasta_line);
            }
            #endregion

            #region read csv file
            // Provide the path to your CSV file
            // Read all lines from the CSV file
            string[] lines = File.ReadAllLines(filePath);

            // Skip the header line (assuming the first line is header)

            ProcessHeaderCSV(lines[0]);
            for (int i = 1; i < lines.Length; i++)
            {
                string _line = Regex.Replace(lines[i], "\"", "");
                // Split the line by comma
                string[] fields = _line.Split(',');

                int index = Array.IndexOf(HeaderLineCSV, "crosslinks_ab");
                if (index == -1) return;
                string crosslinks_ab = fields[index];
                string[] cols_all_xls = Regex.Split(crosslinks_ab, "#");

                index = Array.IndexOf(HeaderLineCSV, "score_ab");
                if (index == -1) return;
                string scores_ab = fields[index];
                string[] cols_all_xls_score = Regex.Split(scores_ab, "#");

                if (cols_all_xls.Length != cols_all_xls_score.Length)
                    return;

                List<Crosslink> crosslinks = new List<Crosslink>();
                for (int count = 0; count < cols_all_xls.Length; count++)
                {
                    string xls = cols_all_xls[count];
                    string[] cols_xl = Regex.Split(xls, "-");
                    Crosslink xl = new Crosslink(cols_xl[0], cols_xl[2], int.Parse(cols_xl[1]), int.Parse(cols_xl[3]), double.Parse(cols_all_xls_score[count]));
                    crosslinks.Add(xl);
                }

                index = Array.IndexOf(HeaderLineCSV, "protein_a");
                if (index == -1) return;
                string ptnA = fields[index];
                if (ptnA.Contains("|"))
                    ptnA = Regex.Split(ptnA, "\\|")[1];


                index = Array.IndexOf(HeaderLineCSV, "protein_b");
                if (index == -1) return;
                string ptnB = fields[index];
                if (ptnB.Contains("|"))
                    ptnB = Regex.Split(ptnB, "\\|")[1];

                FastaItem? fastaA = fasta_items.Where(a => a.accessionNumber.Equals(ptnA)).FirstOrDefault();
                FastaItem? fastaB = fasta_items.Where(a => a.accessionNumber.Equals(ptnB)).FirstOrDefault();

                index = Array.IndexOf(HeaderLineCSV, "length_protein_a");
                if (index == -1) return;
                int lengthA = int.Parse(fields[index]);

                if (fastaA != null && fastaA.length != lengthA)
                {
                    Console.WriteLine($"WARN: Protein {fastaA.accessionNumber} has different length from the informed length of the input csv file. Protein is ignored.");
                    continue;
                }

                index = Array.IndexOf(HeaderLineCSV, "length_protein_b");
                if (index == -1) return;
                int lengthB = int.Parse(fields[index]);

                if (fastaB != null && fastaB.length != lengthB)
                {
                    Console.WriteLine($"WARN: Protein {fastaB.accessionNumber} has different length from the informed length of the input csv file. Protein is ignored.");
                    continue;
                }

                index = Array.IndexOf(HeaderLineCSV, "gene_a");
                if (index == -1) return;
                string geneA = fields[index];

                index = Array.IndexOf(HeaderLineCSV, "gene_b");
                if (index == -1) return;
                string geneB = fields[index];

                Protein proteinA = new();
                proteinA.proteinID = ptnA;
                proteinA.gene = geneA;
                proteinA.fasta = fastaA;
                Protein proteinB = new();
                proteinB.proteinID = ptnB;
                proteinB.gene = geneB;
                proteinB.fasta = fastaB;

                index = Array.IndexOf(HeaderLineCSV, "ppi_score");
                if (index == -1) return;
                double score = double.Parse(fields[index]);

                crosslinks.ForEach(a =>
                {
                    if (a.sourceGene.Equals(proteinA.gene))
                        a.sourceAccessionNumberProtein = ptnA;
                    if (a.targetGene.Equals(proteinB.gene))
                        a.targetAccessionNumberProtein = ptnB;
                });
                crosslinks.RemoveAll(a => a.sourceAccessionNumberProtein == null || a.targetAccessionNumberProtein == null);

                PPI ppi = new(proteinA, proteinB, score, crosslinks);
                ppis.Add(ppi);
            }
            #endregion

            int count_common = 0;
            foreach (var ppi in ppis)
            {
                Management.crosslinks = ppi.crosslinks;
                Management.source_protein_name = ppi.proteinA.proteinID;
                Management.target_protein_name = ppi.proteinB.proteinID;

                Protein ptnA = Management.GetPDBidFromUniprot(ppi.proteinA.proteinID);
                Protein ptnB = Management.GetPDBidFromUniprot(ppi.proteinB.proteinID);

                if (ptnA.proteinID == null || ptnB.proteinID == null) continue;

                List<PDB> common_pdbs = (from pdb_a in ptnA.pdbIDs
                                         from pdb_b in ptnB.pdbIDs
                                         where pdb_a.entry == pdb_b.entry
                                         select pdb_a).ToList();

                if (common_pdbs.Count == 0) continue;

                Management.ProcessPDBFile(common_pdbs[0].entry, ptnA, ptnB, ptnA.gene + "#" + ptnB.gene, false, "", false);
                count_common++;
                Console.WriteLine();
            }
            Console.WriteLine();
        }

        private static void ProcessHeaderCSV(string row)
        {
            row = Regex.Replace(row, "\"", "");
            HeaderLineCSV = Regex.Split(row, ",");
            if (HeaderLineCSV.Length > 1)
            {
                for (int i = 0; i < HeaderLineCSV.Length; i++)
                {
                    HeaderLineCSV[i] = Regex.Replace(HeaderLineCSV[i], "\"", "");
                }

            }
            else
            {
                HeaderLineCSV = Regex.Split(row, ",");
            }
        }
    }
}