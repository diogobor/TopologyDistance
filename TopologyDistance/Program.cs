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
            string mainPath = @"D:\Experimental_data\Cong\";

            //string fileName = $"{mainPath}cytoscapedefaultnode.csv";
            string fileName = $"{mainPath}Scout_DSBSO_cytoscape.csv";
            string fastaFile = $"{mainPath}20220323_CW_DSBSO_Paper_HEK_db.fasta";

            bool useAlphaFold = false;
            int threshold_angstroms = 35;
            string tool = "Scout";
            string crosslinker = "DSBSO";

            string distance_file = $"{mainPath}{tool}_{crosslinker}_distance_output_{threshold_angstroms}A.csv";

            List<PPI> ppis = new();
            if (Management.CleanTmpFiles() == false)
                return;

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
            string[] lines = File.ReadAllLines(fileName);

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


            List<PPI> valid_PPIs = new();

            int ppi_processed = 0;
            int old_progress = 0;
            double ppis_amount = ppis.Count;

            foreach (var ppi in ppis)
            {
                ppi_processed++;
                int new_progress = (int)((double)ppi_processed / (ppis_amount) * 100);
                if (new_progress > old_progress)
                {
                    old_progress = new_progress;
                    int currentLineCursor = Console.CursorTop;
                    Console.SetCursorPosition(0, Console.CursorTop);
                    Console.Write("INFO: PPIs processed: " + old_progress + "%");
                    Console.SetCursorPosition(0, currentLineCursor);
                }

                Management.crosslinks = ppi.crosslinks;
                Management.source_protein_name = ppi.proteinA.proteinID;
                Management.target_protein_name = ppi.proteinB.proteinID;

                if (useAlphaFold)
                {
                    string pdb_source = Management.CreatePDBFile(ppi.proteinA.proteinID, true, false);
                    if (pdb_source.Equals("ERROR"))
                    {
                        Console.WriteLine("ERROR: Unable to create pdb source");
                        continue;
                    }

                    string pdb_target = Management.CreatePDBFile(ppi.proteinB.proteinID, true, true);
                    if (pdb_target.Equals("ERROR"))
                    {
                        Console.WriteLine("ERROR: Unable to create pdb target");
                        continue;
                    }
                    Management.proteinChain_source = "A";
                    Management.HasMoreThanOneChain_proteinTarget = false;

                    bool valid = Management.ProcessPDBorCIFfileWithSpecificChain(ppi.proteinA, ppi.proteinB, "B", pdb_source, pdb_target, true);
                    if (valid)
                        valid_PPIs.Add(ppi);

                }
                else
                {
                    Protein ptnA = Management.GetPDBidFromUniprot(ppi.proteinA.proteinID);
                    Protein ptnB = Management.GetPDBidFromUniprot(ppi.proteinB.proteinID);

                    if (ptnA.proteinID == null || ptnB.proteinID == null) continue;

                    List<PDB> common_pdbs = (from pdb_a in ptnA.pdbIDs
                                             from pdb_b in ptnB.pdbIDs
                                             where pdb_a.entry == pdb_b.entry
                                             select pdb_a).ToList();

                    if (common_pdbs.Count == 0) continue;

                    #region check cross-links positions and protein length
                    List<Crosslink> crosslinks_ptn = ppi.crosslinks;

                    string[] cols = Regex.Split(ptnA.pdbIDs[0].positions, "-");
                    int min_pos_pdb_a = int.Parse(cols[0]);
                    int max_pos_pdb_a = -1;
                    if (cols[1].Contains(","))
                        max_pos_pdb_a = int.Parse(Regex.Split(cols[1], ",")[0]);
                    else
                        max_pos_pdb_a = int.Parse(cols[1]);

                    cols = Regex.Split(ptnB.pdbIDs[0].positions, "-");
                    int min_pos_pdb_b = int.Parse(cols[0]);
                    int max_pos_pdb_b = -1;
                    if (cols[1].Contains(","))
                        max_pos_pdb_b = int.Parse(Regex.Split(cols[1], ",")[0]);
                    else
                        max_pos_pdb_b = int.Parse(cols[1]);

                    crosslinks_ptn.RemoveAll(a => a.sourceAccessionNumberProtein.Equals(ppi.proteinA.proteinID) &&
                    (a.sourcePosition + min_pos_pdb_a) > max_pos_pdb_a);

                    crosslinks_ptn.RemoveAll(a => a.targetAccessionNumberProtein.Equals(ppi.proteinA.proteinID) &&
                    (a.targetPosition + min_pos_pdb_a) > max_pos_pdb_a);

                    crosslinks_ptn.RemoveAll(a => a.sourceAccessionNumberProtein.Equals(ppi.proteinB.proteinID) &&
                    (a.sourcePosition + min_pos_pdb_b) > max_pos_pdb_b);

                    crosslinks_ptn.RemoveAll(a => a.targetAccessionNumberProtein.Equals(ppi.proteinB.proteinID) &&
                    (a.targetPosition + min_pos_pdb_b) > max_pos_pdb_b);

                    if (crosslinks_ptn.Count == 0) continue;

                    #endregion


                    bool valid = Management.ProcessPDBFile(common_pdbs[0].entry, ptnA, ptnB, ptnA.gene + "#" + ptnB.gene, false, "", false, false);
                    if (valid)
                        valid_PPIs.Add(ppi);
                }
                
                Console.WriteLine();
            }

            Console.WriteLine("Done!");
            Console.WriteLine($"{valid_PPIs.Count} from {ppis.Count} have been computed.");

            int lower_threshold = 0;
            int upper_threshold = 0;
            StringBuilder sb_distances = new();
            foreach (var _ppi in valid_PPIs)
            {
                lower_threshold += _ppi.crosslinks.Count(a => a.distance <= threshold_angstroms);
                upper_threshold += _ppi.crosslinks.Count(a => a.distance > threshold_angstroms);
                sb_distances.AppendLine(String.Join('\n', _ppi.crosslinks.Select(a => a.distance).ToList()));
            }

            Console.WriteLine($"Total XLs: {valid_PPIs.Sum(a => a.crosslinks.Count)}");
            Console.WriteLine($"Lower than {threshold_angstroms}A: {lower_threshold}\nUpper than {threshold_angstroms}A: {upper_threshold}");

            File.WriteAllText(distance_file, sb_distances.ToString());
            Management.CleanTmpFiles();
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