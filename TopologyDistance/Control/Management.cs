using System;
using System.Collections;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Net;
using System.Runtime.CompilerServices;
using System.Text;
using System.Text.RegularExpressions;
using System.Threading.Tasks;
using System.Xml;
using TopologyDistance.Model;

namespace TopologyDistance.Control
{
    public static class Management
    {
        private static int REQUEST_TIMEOUT = 10000;
        private static int proteinOffsetInPDBTarget { get; set; } = -1;
        private static int proteinOffsetInPDBSource { get; set; } = -1;
        private static string proteinSequenceFromPDBFile_proteinTarget;
        private static bool HasMoreThanOneChain_proteinSource;
        private static string proteinSequenceFromPDBFile_proteinSource;
        private static string pdbFile;
        private static Aligner align = new();

        public static bool HasMoreThanOneChain_proteinTarget;
        public static string proteinChain_source;
        public static string source_protein_name;
        public static string target_protein_name;
        public static List<Crosslink> crosslinks { get; set; }

        public static Protein GetPDBidFromUniprot(string proteinID)
        {
            try
            {
                string responseString = GetProteinInfoFromUniprot(proteinID);
                if (!string.IsNullOrWhiteSpace(responseString))
                {
                    if (responseString.StartsWith("<!DOCTYPE html PUBLIC"))
                        return new Protein();
                    XmlDocument doc = new XmlDocument();
                    doc.LoadXml(responseString);
                    if (doc == null)
                        return new Protein();

                    XmlNodeList xmlnodes = doc.GetElementsByTagName("error");
                    if (xmlnodes.Count > 0)
                    {
                        throw new Exception("ERROR: " + xmlnodes[0].InnerText);
                    }

                    //Console.WriteLine("INFO: Getting protein description...");
                    xmlnodes = doc.GetElementsByTagName("recommendedName");
                    if (xmlnodes.Count == 0)
                    {
                        xmlnodes = doc.GetElementsByTagName("submittedName");
                    }

                    XmlNodeList nodes = xmlnodes[0].ChildNodes;
                    string fullName = "";
                    foreach (XmlNode node in nodes)
                    {
                        if (node is XmlElement)
                        {
                            if (node.Name.Equals("fullName"))
                            {
                                fullName = node.FirstChild.InnerText;
                                break;
                            }
                        }
                    }

                    //Console.WriteLine("INFO: Getting gene name...");

                    xmlnodes = doc.GetElementsByTagName("gene");
                    nodes = xmlnodes[0].ChildNodes;
                    string geneName = "";
                    foreach (XmlNode node in nodes)
                    {
                        if (node is XmlElement)
                        {
                            if (node.Name.Equals("name"))
                            {
                                geneName = node.FirstChild.InnerText;
                                break;
                            }
                        }
                    }

                    //Console.WriteLine("INFO: Getting PDB IDs...");

                    xmlnodes = doc.GetElementsByTagName("dbReference");
                    bool containsPDBtags = false;
                    List<PDB> pdbs = new List<PDB>();
                    foreach (XmlNode nNode in xmlnodes)
                    {
                        string pdbType = nNode.Attributes[0].Value;
                        if (pdbType.Equals("PDB"))
                        {
                            if (nNode.HasChildNodes)
                            {
                                string entry = nNode.Attributes[1].Value;
                                string resolution = "0.00";
                                string[] chain_positions = null;
                                if (nNode.ChildNodes[1].Attributes[0].Value.Equals("resolution"))
                                {
                                    resolution = nNode.ChildNodes[1].Attributes[1].Value;
                                    chain_positions = nNode.ChildNodes[2].Attributes[1].Value.Split('=');
                                }
                                else if (nNode.ChildNodes[1].Attributes[0].Value.Equals("chains"))
                                {
                                    chain_positions = nNode.ChildNodes[1].Attributes[1].Value.Split('=');
                                }
                                string chain = chain_positions[0];
                                string positions = chain_positions[1];
                                PDB pdb = new PDB(entry, resolution, chain, positions);
                                pdbs.Add(pdb);
                                containsPDBtags = true;
                            }
                        }
                        else if (pdbType.Equals("SMR") && !containsPDBtags)
                        {
                            string entry = nNode.Attributes[1].Value;
                            PDB pdb = new PDB(entry, "SMR", "", "");
                            pdbs.Add(pdb);
                        }
                    }

                    //Console.WriteLine("INFO: Getting protein sequence...");

                    xmlnodes = doc.GetElementsByTagName("sequence");
                    string ptnSequence = "";
                    string checksum = "";
                    foreach (XmlNode node in xmlnodes)
                    {
                        if (node is XmlElement)
                        {
                            if (node.Attributes != null && node.Attributes.Count > 2 && node.Attributes[2].Name.Equals("checksum"))
                            {
                                checksum = node.Attributes[2].Value;
                                ptnSequence = node.FirstChild.InnerText;
                                break;
                            }
                        }
                    }

                    Protein ptn = new Protein(proteinID, fullName, geneName, ptnSequence, checksum, pdbs, null);
                    return ptn;
                }
                else
                {
                    return new Protein();
                }
            }
            catch (Exception e)
            {
                Console.WriteLine(e.Message);
                return new Protein();
            }
        }
        public static string CreatePDBFile(string pdbID, bool isAlphaFold, bool modifyChain)
        {
            string finalStr = "";
            string[] returnFile = null;
            FileInfo f = null;
            try
            {
                if (isAlphaFold)
                {
                    // Download file from AlphaFold server
                    returnFile = GetPDBfileFromAlphaFoldServer(pdbID, modifyChain);
                }
                else
                {
                    if (pdbID.StartsWith("https://swissmodel.expasy.org/repository/"))
                    {
                        // Download file from SwissModel server
                        returnFile = GetPDBfileFromSwissModelServer(pdbID);
                    }
                    else
                    {
                        // Retrieving file from RCSB server
                        // [PDB or CIF, file content]
                        returnFile = GetPDBorCIFfileFromServer(pdbID);
                    }
                }

                if (String.IsNullOrEmpty(returnFile[0]))
                {
                    Console.WriteLine("ERROR: Unable to retrieve protein structure");
                    return "";
                }
                // write this script to tmp file and return path
                if (returnFile[0].Equals("PDB"))
                    f = GetTmpFile(pdbID, "pdb");
                else
                    f = GetTmpFile(pdbID, "cif");

                if (f == null)
                {
                    Console.WriteLine("ERROR: Could not create temp file with label: " + pdbID);
                    return "";
                }
                using (StreamWriter bw = new StreamWriter(f.FullName))
                {
                    finalStr = returnFile[1];
                    if (string.IsNullOrWhiteSpace(finalStr))
                        Console.WriteLine("ERROR: Error retrieving the PDB file for protein " + pdbID + ".");
                    bw.Write(finalStr);
                }
            }
            catch (IOException ex)
            {
                Console.WriteLine("ERROR: Problems while writing script to " + f.FullName);
                finalStr = "";
            }
            if (string.IsNullOrWhiteSpace(finalStr))
                return "ERROR";
            return f != null ? f.FullName : "ERROR";
        }
        public static string GetChainFromPDBFasta(Protein ptn, string pdbID, string fileName)
        {
            List<FastaItem> fastaList = GetProteinSequenceFromPDBServer(pdbID);
            if (fastaList.Count == 0)
            {
                fastaList = GetProteinSequencesFromPDBFile(fileName, ptn);
            }
            string chain = "";
            if (fastaList.Count > 0)
            {
                foreach (FastaItem fasta in fastaList)
                {
                    int offset = fasta.sequence.IndexOf(ptn.fasta.sequence);
                    if (offset == -1)
                        offset = ptn.fasta.sequence.IndexOf(fasta.sequence);
                    if (offset == -1)
                    {
                        try
                        {
                            (int MaxScorePos, string closestPeptideInProteinDB) closestPeptideinProteinInfo = align.GetClosestPeptideInASequence(fasta.sequence.ToCharArray(), ptn.fasta.sequence.ToCharArray());
                            string closestPept = (string)closestPeptideinProteinInfo.closestPeptideInProteinDB;
                            (int MaxScorePos, string closestPeptideInProteinDB) closestPeptInfo;
                            if (fasta.sequence.Length > closestPept.Length)
                                closestPeptInfo = align.GetClosestPeptideInASequence(closestPept.ToCharArray(), fasta.sequence.ToCharArray());
                            else
                                closestPeptInfo = align.GetClosestPeptideInASequence(fasta.sequence.ToCharArray(), closestPept.ToCharArray());
                            int indexClosestPet = closestPeptInfo.MaxScorePos;
                            int countSameAA = 0;
                            int limit = Math.Min(closestPept.Length, fasta.sequence.Length - indexClosestPet);
                            for (int i = indexClosestPet; i < limit; i++)
                            {
                                if (closestPept[i] == fasta.sequence[i])
                                    countSameAA++;
                            }
                            if (((double)countSameAA / (double)limit) > 0.4)
                                offset = 0;
                        }
                        catch (Exception e)
                        {
                            continue;
                        }
                    }
                    if (offset != -1)
                    {
                        // Example: >3J7Y_8|Chain J|uL11|Homo sapiens (9606)
                        string[] cols = fasta.accessionNumber.Split('|');
                        string[] chainCols = cols[1].Split("Chain ");
                        if (chainCols.Length == 1)
                        {
                            chainCols = cols[1].Split("Chains ");
                            chainCols = chainCols[1].Split(',');
                            chain = chainCols[0][0].ToString();
                        }
                        else
                        {
                            chain = chainCols[1][0].ToString();
                        }
                        if (proteinOffsetInPDBSource == -1 && fasta.offset != -1)
                            proteinOffsetInPDBSource = fasta.offset;
                        break;
                    }
                }
            }
            return chain;
        }
        public static string[] GetProteinSequenceAndChainFromPDBFile(string fileName, Protein ptn)
        {
            Dictionary<string, int> ResiduesDict = CreateResiduesDict();
            StringBuilder sbSequence = new StringBuilder();
            string protein_chain = "";
            StringBuilder sbProteinChains = new StringBuilder();
            sbProteinChains.Append("CHAINS:");
            bool hasMoreThanOneChain = false;
            int countChains = 0;
            try
            {
                using (StreamReader parserFile = new StreamReader(fileName))
                {
                    string line;
                    int lastInsertedResidue = 0;
                    int threshold = 10; // qtd aminoacids
                    int countAA = 0;
                    bool isCompleteFullName = false;
                    StringBuilder sbProteinFullName = new StringBuilder();
                    while ((line = parserFile.ReadLine()) != null)
                    {
                        if (string.IsNullOrEmpty(line))
                            continue;
                        if (!line.StartsWith("COMPND") && !line.StartsWith("ATOM"))
                            continue;
                        if (line.StartsWith("COMPND"))
                        {
                            string[] cols = line.Split(new[] { ' ' }, StringSplitOptions.RemoveEmptyEntries);
                            if (double.TryParse(cols[1], out double isNumeric))
                            {
                                // Get protein full name
                                if (cols[2].Equals("MOLECULE:"))
                                {
                                    sbProteinFullName.Clear();
                                    for (int i = 3; i < cols.Length; i++)
                                    {
                                        sbProteinFullName.Append(cols[i] + " ");
                                    }
                                    if (cols[^1].EndsWith(";"))
                                    {
                                        isCompleteFullName = true;
                                        continue;
                                    }
                                }
                                else if (!isCompleteFullName)
                                {
                                    sbProteinFullName.Clear();
                                    for (int i = 2; i < cols.Length; i++)
                                    {
                                        sbProteinFullName.Append(cols[i] + " ");
                                    }
                                    if (cols[^1].EndsWith(";"))
                                    {
                                        isCompleteFullName = true;
                                        continue;
                                    }
                                }
                                if (isCompleteFullName)
                                {
                                    if (cols[2].Equals("CHAIN:"))
                                    {
                                        countChains++;
                                        string pdbProteinFullName_moleculeField = sbProteinFullName.ToString()
                                            .Replace(';', ' ').Trim();
                                        string[] fullNames = pdbProteinFullName_moleculeField.Split(',');
                                        if (fullNames.Length == 1)
                                        {
                                            if (ptn.fullName.ToLower()
                                                .Equals(pdbProteinFullName_moleculeField.ToLower()))
                                            {
                                                protein_chain = cols[3].ToString().Replace(';', ' ').Replace(',', ' ')
                                                    .Trim();
                                                continue;
                                            }
                                        }
                                        else
                                        {
                                            // There is only one molecule, but with more than one chain for this ptn full name
                                            hasMoreThanOneChain = true;
                                            foreach (string fullName in fullNames)
                                            {
                                                if (ptn.fullName.ToLower().Equals(fullName.ToLower().Trim()))
                                                {
                                                    protein_chain = cols[3].ToString().Replace(';', ' ').Trim()
                                                        .Split(",")[0].Trim();
                                                    break;
                                                }
                                            }
                                        }
                                        // Check if 'CHAIN' field has more than 1 chain: CHAIN: A,B
                                        if (cols.Length - 3 > 1)
                                            hasMoreThanOneChain = true;
                                        for (int i = 3; i < cols.Length; i++)
                                        {
                                            sbProteinChains.Append(
                                                cols[i].ToString().Replace(',', ' ').Replace(';', ' ').Trim() + "#");
                                        }
                                    }
                                }
                            }
                            else
                            {
                                continue;
                            }
                        }
                        else // It starts with 'ATOM'
                        {
                            string[] cols = line.Split(new[] { ' ' }, StringSplitOptions.RemoveEmptyEntries);
                            if (!cols[4].Equals(protein_chain))
                                continue;
                            string pdbResidue = cols[3]; // Residue -> three characters
                            int newResidue = -1;
                            ResiduesDict.TryGetValue(pdbResidue, out newResidue);

                            if (newResidue != lastInsertedResidue)
                            {
                                sbSequence.Append((char)newResidue);
                                countAA++;
                                if (countAA > threshold)
                                {
                                    break;
                                }
                            }
                            lastInsertedResidue = newResidue;
                            if (proteinOffsetInPDBSource == -1)
                            {
                                proteinOffsetInPDBSource = int.Parse(cols[5]);
                            }
                        }
                    }
                }
            }
            catch (Exception e)
            {
                Console.WriteLine("ERROR: Problems while reading PDB file: " + fileName);
            }

            if (string.IsNullOrWhiteSpace(protein_chain))
                protein_chain = sbProteinChains.ToString();

            if (countChains > 1 || hasMoreThanOneChain)
            {
                return new string[] { sbSequence.ToString(), protein_chain, "true" };
            }
            else
            {
                return new string[] { sbSequence.ToString(), protein_chain, "false" };
            }
        }
        public static string[] GetChainFromCIFFile(string fileName)
        {
            StringBuilder sbSequence = new StringBuilder();
            string protein_chain = "";
            HashSet<string> chainsSet = new HashSet<string>();
            StringBuilder sbProteinChains = new StringBuilder();
            sbProteinChains.Append("CHAINS:");
            int countChains = 0;

            try
            {
                using (StreamReader parserFile = new StreamReader(fileName))
                {
                    string line;
                    int proteinOffsetInPDBSource = -1;
                    while ((line = parserFile.ReadLine()) != null)
                    {
                        if (string.IsNullOrEmpty(line))
                            continue;
                        if (!line.StartsWith("ATOM"))
                            continue;
                        string[] cols = line.Split(new[] { ' ' }, StringSplitOptions.RemoveEmptyEntries);
                        if (cols[5].Length == 3) // It means that the ATOM is a residue, not a gene
                            chainsSet.Add(cols[6]);
                    }
                }
            }
            catch (Exception e)
            {
                Console.WriteLine("ERROR: Problems while reading PDB file: " + fileName);
            }

            countChains = sbProteinChains.Length;
            sbProteinChains.Append(string.Join("#", chainsSet));
            protein_chain = sbProteinChains.ToString();

            if (countChains > 1)
            {
                return new string[] { sbSequence.ToString(), protein_chain, "true" };
            }
            else
            {
                return new string[] { sbSequence.ToString(), protein_chain, "false" };
            }
        }

        //python
        private static List<string> OutputLog { get; set; } = new();
        private static List<string> ErrorLog { get; set; } = new();
        private static EventHandler OnExit { get; set; }

        public static bool ProcessPDBFile(
            string pdbID,
            Protein ptnSource,
            Protein ptnTarget,
            string nodeName,
            bool processTarget,
            string proteinChain_proteinSource,
            bool useCustomizedPDBFile,
            bool useAlphaFold)
        {
            if (processTarget)
            {
                //Console.WriteLine("INFO: Chain of protein source: " + pdbID);
                proteinChain_source = proteinChain_proteinSource;
                HasMoreThanOneChain_proteinTarget = false;
                bool foundChain = true;
                string proteinChain_proteinTarget = GetChainFromPDBFasta(ptnTarget, pdbID, pdbFile);
                if (string.IsNullOrWhiteSpace(proteinChain_proteinTarget))
                    foundChain = false;
                if (!foundChain)
                {
                    string[] returnPDB_proteinTarget = null;
                    if (pdbFile.EndsWith("pdb"))
                    {
                        //Console.WriteLine("INFO: Getting protein sequence and chain of protein source from PDB file...");
                        returnPDB_proteinTarget = GetProteinSequenceAndChainFromPDBFile(pdbFile, ptnTarget);
                    }
                    else
                    {
                        //Console.WriteLine("INFO: Getting all chains of protein source from CIF file...");
                        returnPDB_proteinTarget = GetChainFromCIFFile(pdbFile);
                    }
                    proteinSequenceFromPDBFile_proteinTarget = returnPDB_proteinTarget[0];
                    HasMoreThanOneChain_proteinTarget = returnPDB_proteinTarget[2] == "true";
                    proteinChain_proteinTarget = returnPDB_proteinTarget[1];
                    if (proteinChain_proteinTarget.StartsWith("CHAINS:"))
                    {
                        //Console.WriteLine("WARN: No chain matched with protein target description.");
                        string[] protein_chains = returnPDB_proteinTarget[1].Replace("CHAINS:", "").Split('#');
                        List<PDB> PDBchains = new List<PDB>();
                        foreach (string chainStr in protein_chains)
                        {
                            PDBchains.Add(new PDB("", "", chainStr, ""));
                        }
                        if (PDBchains.Count > 1)
                        {
                            //Console.WriteLine("INFO: There is more than one chain. The first one was selected...");
                            // Open a window to select only one chain
                            return ProcessPDBorCIFfileWithSpecificChain(ptnSource, ptnTarget, PDBchains[0].chain, pdbFile, pdbFile, useAlphaFold);
                        }
                    }
                    else
                    {
                        return ProcessPDBorCIFfileWithSpecificChain(ptnSource, ptnTarget, returnPDB_proteinTarget[1], pdbFile, pdbFile, useAlphaFold);
                    }
                }
                else
                {
                    return ProcessPDBorCIFfileWithSpecificChain(ptnSource, ptnTarget, proteinChain_proteinTarget, pdbFile, pdbFile, useAlphaFold);
                }
            }
            else
            {
                //Console.WriteLine("INFO: Creating tmp PDB file...");
                if (!useCustomizedPDBFile)
                {
                    pdbFile = CreatePDBFile(pdbID, useAlphaFold, false);
                }
                if (string.IsNullOrWhiteSpace(pdbFile) || pdbFile.Equals("ERROR"))
                {
                    Console.WriteLine("ERROR: Error creating PDB file.");
                    return false;
                }
                HasMoreThanOneChain_proteinSource = false;
                bool foundChain = true;
                proteinChain_proteinSource = GetChainFromPDBFasta(ptnSource, pdbID, pdbFile);
                if (string.IsNullOrWhiteSpace(proteinChain_proteinSource))
                    foundChain = false;
                if (!foundChain)
                {
                    string[] returnPDB_proteinSource = null;
                    if (!string.IsNullOrWhiteSpace(pdbFile) && pdbFile.EndsWith("pdb"))
                    {
                        //Console.WriteLine("INFO: Getting protein sequence and chain of protein source from PDB file...");
                        returnPDB_proteinSource = GetProteinSequenceAndChainFromPDBFile(pdbFile, ptnSource);
                    }
                    else
                    {
                        //Console.WriteLine("INFO: Getting all chains of protein source from CIF file...");
                        returnPDB_proteinSource = GetChainFromCIFFile(pdbFile);
                    }
                    proteinSequenceFromPDBFile_proteinSource = returnPDB_proteinSource[0];
                    HasMoreThanOneChain_proteinSource = returnPDB_proteinSource[2] == "true";
                    proteinChain_proteinSource = returnPDB_proteinSource[1];
                    if (proteinChain_proteinSource.StartsWith("CHAINS:"))
                    {
                        //Console.WriteLine("WARN: No chain matched with protein source description.");
                        string[] protein_chains = returnPDB_proteinSource[1].Replace("CHAINS:", "").Split('#');
                        List<PDB> PDBchains = new List<PDB>();
                        foreach (string chainStr in protein_chains)
                        {
                            PDBchains.Add(new PDB("", "", chainStr, ""));
                        }
                        if (PDBchains.Count > 1)
                        {
                            //Console.WriteLine("INFO: There is more than one chain. Select first one was selected...");
                            // Open a window to select only one chain
                            ProcessPDBFile(pdbID, ptnSource, ptnTarget, nodeName, true, PDBchains[0].chain, useCustomizedPDBFile, useAlphaFold);
                        }
                    }
                    else
                    {
                        ProcessPDBFile(pdbID, ptnSource, ptnTarget, nodeName, true, proteinChain_proteinSource, useCustomizedPDBFile, useAlphaFold);
                    }
                }
                else
                {
                    ProcessPDBFile(pdbID, ptnSource, ptnTarget, nodeName, true, proteinChain_proteinSource, useCustomizedPDBFile, useAlphaFold);
                }
            }
            return true;
        }

        public static string GetProteinSequenceFromPDBFileWithSpecificChain(
            string fileName,
            string protein_chain,
            bool isProteinSource)
        {
            var ResiduesDict = CreateResiduesDict();
            var sbSequence = new StringBuilder();
            int proteinOffsetInPDB = -1;

            try
            {
                using (StreamReader parserFile = new StreamReader(fileName))
                {
                    string line;
                    string lastOffset = "-1";
                    int threshold = 15; // Number of amino acids
                    int countAA = 0;
                    bool getSequence = true;

                    while ((line = parserFile.ReadLine()) != null)
                    {
                        if (!string.IsNullOrEmpty(line))
                        {
                            if (!line.StartsWith("ATOM"))
                                continue;
                            // It starts with 'ATOM'
                            var cols = line.Split(new char[] { ' ', '\t' }, StringSplitOptions.RemoveEmptyEntries);
                            if (cols[4] != protein_chain)
                                continue;
                            var pdbResidue = cols[3]; // Residue -> three characters

                            int newResidue = -1;
                            ResiduesDict.TryGetValue(pdbResidue, out newResidue);

                            var currentOffset = cols[5];
                            if (currentOffset != lastOffset)
                            {
                                var _byte = new byte[1];
                                _byte[0] = (byte)newResidue;
                                var _string = Encoding.UTF8.GetString(_byte);
                                sbSequence.Append(_string);
                                countAA++;
                                if (countAA > threshold)
                                {
                                    break;
                                }
                            }
                            lastOffset = currentOffset;
                            if (proteinOffsetInPDB == -1)
                            {
                                proteinOffsetInPDB = int.Parse(cols[5]);
                            }
                        }
                    }
                }
            }
            catch (Exception e)
            {
                Console.WriteLine("ERROR: Problems while reading PDB file: " + fileName);
            }

            if (isProteinSource)
                proteinOffsetInPDBSource = proteinOffsetInPDB;
            else
                proteinOffsetInPDBTarget = proteinOffsetInPDB;

            return sbSequence.ToString();
        }

        public static string GetProteinSequenceFromCIFFileWithSpecificChain(
            string fileName,
            string protein_chain,
            bool isProteinSource)
        {
            var ResiduesDict = CreateResiduesDict();
            var sbSequence = new StringBuilder();
            int proteinOffsetInPDB = -1;

            try
            {
                using (StreamReader parserFile = new StreamReader(fileName))
                {
                    string line;
                    string lastOffset = "-1";
                    int threshold = 15; // Number of amino acids
                    int countAA = 0;
                    bool getSequence = true;

                    while ((line = parserFile.ReadLine()) != null)
                    {
                        if (!string.IsNullOrEmpty(line))
                        {
                            if (!line.StartsWith("ATOM"))
                                continue;
                            // It starts with 'ATOM'
                            var cols = line.Split(new char[] { ' ', '\t' }, StringSplitOptions.RemoveEmptyEntries);
                            if (cols[6] != protein_chain)
                                continue;
                            if (cols[5].Length != 3)
                                continue; // It means that the ATOM is not a residue, it's a gene
                            var pdbResidue = cols[5]; // Residue -> three characters

                            int newResidue = -1;
                            ResiduesDict.TryGetValue(pdbResidue, out newResidue);

                            var currentOffset = cols[7];
                            if (currentOffset != lastOffset)
                            {
                                var _byte = new byte[1];
                                _byte[0] = (byte)newResidue;
                                var _string = Encoding.UTF8.GetString(_byte);
                                sbSequence.Append(_string);
                                countAA++;
                                if (countAA > threshold)
                                {
                                    break;
                                }
                            }
                            lastOffset = currentOffset;
                            if (proteinOffsetInPDB == -1)
                            {
                                proteinOffsetInPDB = int.Parse(cols[7]);
                            }
                        }
                    }
                }
            }
            catch (Exception e)
            {
                Console.WriteLine("ERROR: Problems while reading PDB file: " + fileName);
            }

            proteinOffsetInPDBSource = -1;
            proteinOffsetInPDBTarget = -1;
            if (isProteinSource)
                proteinOffsetInPDBSource = proteinOffsetInPDB;
            else
                proteinOffsetInPDBTarget = proteinOffsetInPDB;

            return sbSequence.ToString();
        }

        public static bool ProcessPDBorCIFfileWithSpecificChain(
            Protein ptnSource,
            Protein ptnTarget,
            string proteinChain_target,
            string pdbFile_source,
            string pdbFile_target,
            bool useAlphaFold)
        {
            //Console.WriteLine("INFO: Chain of protein target: " + proteinChain_target);
            //Console.WriteLine("INFO: Getting sequence of protein source: " + source_protein_name);

            string proteinSequence_source_FromPDBFile = "";
            if (pdbFile_source.EndsWith("pdb"))
                proteinSequence_source_FromPDBFile = GetProteinSequenceFromPDBFileWithSpecificChain(pdbFile_source, proteinChain_source, true);
            else
                proteinSequence_source_FromPDBFile = GetProteinSequenceFromCIFFileWithSpecificChain(pdbFile_source, proteinChain_source, true);

            if (string.IsNullOrEmpty(proteinSequence_source_FromPDBFile))
            {
                Console.WriteLine("ERROR: No sequence has been found in pdb/cif file for: " + source_protein_name);
                return false;
            }

            //Console.WriteLine("INFO: Getting sequence of protein target: " + target_protein_name);

            string proteinSequence_target_FromPDBFile = "";
            if (pdbFile_target.EndsWith("pdb"))
                proteinSequence_target_FromPDBFile = GetProteinSequenceFromPDBFileWithSpecificChain(pdbFile_target, proteinChain_target, false);
            else
                proteinSequence_target_FromPDBFile = GetProteinSequenceFromCIFFileWithSpecificChain(pdbFile_target, proteinChain_target, false);

            if (string.IsNullOrEmpty(proteinSequence_target_FromPDBFile))
            {
                Console.WriteLine("ERROR: No sequence has been found in pdb/cif file for: " + target_protein_name);
                return false;
            }

            // Filter cross-links to obtain only links that belong to source and target nodes
            crosslinks.RemoveAll(o => !(o.sourceAccessionNumberProtein.Equals(source_protein_name) && o.targetAccessionNumberProtein.Equals(target_protein_name)
                                    || o.sourceAccessionNumberProtein.Equals(target_protein_name) && o.targetAccessionNumberProtein.Equals(source_protein_name)));

            string tmpPyMOLScriptFile = CreatePyMOLScriptFile(
                ptnSource,
                ptnTarget,
                crosslinks,
                pdbFile_source,
                pdbFile_target,
                proteinSequence_source_FromPDBFile,
                proteinSequence_target_FromPDBFile,
                HasMoreThanOneChain_proteinSource,
                HasMoreThanOneChain_proteinTarget,
                proteinChain_source,
                proteinChain_target,
                source_protein_name,
                target_protein_name,
                useAlphaFold);

            if (tmpPyMOLScriptFile.Equals("ERROR"))
            {
                Console.WriteLine("ERROR: Error creating PyMOL script file.");
                throw new Exception("Error creating PyMOL script file.");
            }

            ExecutePyMOL(tmpPyMOLScriptFile);

            if ((OutputLog.FindIndex(a => a.StartsWith("ERROR: There is no")) != -1) ||
                ErrorLog.Count > 0) return false;
            string outputfile = OutputLog.Find(a => a.StartsWith("[[")).Trim();

            //Manage outputfile distances
            ManageXLDistances(outputfile);

            return true;
        }

        public static bool CleanTmpFiles()
        {
            try
            {
                DirectoryInfo dr = new DirectoryInfo(Path.GetTempPath() + "\\cytoTmpScripts");
                if (dr.Exists)
                {
                    string directoryPath = dr.FullName;

                    // Delete *.py files
                    DeleteFiles(directoryPath, "*.py");
                    // Delete *.pml files
                    DeleteFiles(directoryPath, "*.pml");
                    // Delete *.pdb files
                    DeleteFiles(directoryPath, "*.pdb");
                    // Delete *.cif files
                    DeleteFiles(directoryPath, "*.cif");
                }
                return true;
            }
            catch (Exception)
            {
                Console.WriteLine("ERROR: Unable to delete files.");
                return false;
            }
        }
        internal static void DeleteFiles(string directoryPath, string searchPattern)
        {
            try
            {
                // Get files matching the search pattern
                string[] files = Directory.GetFiles(directoryPath, searchPattern);
                // Delete each file
                foreach (string file in files)
                    File.Delete(file);
            }
            catch (Exception ex)
            {
                Console.WriteLine($"An error occurred: {ex.Message}");
            }
        }

        internal static void ExecutePyMOL(string script_file)
        {
            OutputLog.Clear();
            ErrorLog.Clear();

            Process proc = new Process();
            proc.StartInfo.FileName = @"C:\ProgramData\pymol\python.exe";
            proc.StartInfo.RedirectStandardOutput = true;
            proc.StartInfo.RedirectStandardError = true;
            proc.StartInfo.UseShellExecute = false;
            proc.StartInfo.CreateNoWindow = true;
            proc.StartInfo.Arguments = script_file.Trim();
            proc.StartInfo.WorkingDirectory = AppDomain.CurrentDomain.BaseDirectory;

            proc.OutputDataReceived += UpdateOutputLog;
            proc.ErrorDataReceived += UpdateErrorLog;

            proc.OutputDataReceived += PrintToConsoleEvent;
            proc.ErrorDataReceived += PrintToConsoleEvent;
            proc.Exited += OnExit;

            proc.Start();
            proc.BeginOutputReadLine();
            proc.BeginErrorReadLine();
            proc.WaitForExit();
        }

        internal static void UpdateOutputLog(object sender, DataReceivedEventArgs e)
        {
            if (e.Data != null)
            {
                OutputLog.Add(e.Data);
            }
        }

        internal static void UpdateErrorLog(object sender, DataReceivedEventArgs e)
        {
            if (e.Data != null)
            {
                ErrorLog.Add(e.Data);
            }
        }

        internal static void PrintToConsoleEvent(object sender, DataReceivedEventArgs e)
        {
            if (e.Data != null)
                Console.WriteLine(e.Data);
        }

        internal static void ManageXLDistances(string distances)
        {
            // Read the text from the file
            //string text = File.ReadAllText(distances);

            // Parse the text into a list of tuples
            List<(int, double)> dataList = ParseData(distances);

            if (crosslinks.Count != dataList.Count)
            {
                Console.WriteLine("ERROR: Crosslink and dataList size are different.");
                return;
            }
            for (int i = 0; i < crosslinks.Count; i++)
            {
                Crosslink current_xl = crosslinks[i];
                current_xl.distance = dataList[i].Item2;
            }

            //File.Delete(distances);
        }
        internal static List<(int, double)> ParseData(string text)
        {
            List<(int, double)> dataList = new();

            var matches = Regex.Matches(text, "(\\[xl)(\\d+)(, )(\\d+.\\d+)(\\])");
            if (matches != null)
            {
                foreach (Match match in matches)
                {
                    GroupCollection gc = match.Groups;
                    dataList.Add((Convert.ToInt32(gc[2].Value), Convert.ToDouble(gc[4].Value)));
                }
            }

            return dataList;
        }
        internal static string CreatePyMOLScriptFile(
            Protein ptnSource,
            Protein ptnTarget,
            List<Crosslink> crossLinks,
            string pdbFile_source,
            string pdbFile_target,
            string proteinSequence_source_FromPDBFile,
            string proteinSequence_target_FromPDBFile,
            bool HasMoreThanOneChain_proteinSource,
            bool HasMoreThanOneChain_proteinTarget,
            string proteinChain_source,
            string proteinChain_target,
            string nodeName_source,
            string nodeName_target,
            bool useAlphaFold)
        {
            string finalStr = "";
            FileInfo filePath = null;
            try
            {
                // Write this script to tmp file and return path
                filePath = GetTmpFile(ptnSource.proteinID + "_" + ptnTarget.proteinID, "py");
                //filePath = GetTmpFile(ptnSource.proteinID + "_" + ptnTarget.proteinID, "pml");

                if (filePath == null)
                {
                    Console.WriteLine("ERROR: Could not create temp file with label: " + ptnSource.proteinID + "_" + ptnTarget.proteinID);
                    return "";
                }


                if (string.IsNullOrEmpty(proteinSequence_source_FromPDBFile) || string.IsNullOrEmpty(proteinSequence_target_FromPDBFile))
                {
                    Console.WriteLine("ERROR: No sequence has been found in PDB file.");
                    return "ERROR";
                }
                using (StreamWriter bw = new StreamWriter(filePath.FullName))
                {
                    finalStr = CreatePythonScript(
                    //finalStr = CreatePyMOLScript(
                        ptnSource,
                        ptnTarget,
                        crossLinks,
                        pdbFile_source,
                        pdbFile_target,
                        HasMoreThanOneChain_proteinSource,
                        HasMoreThanOneChain_proteinTarget,
                        proteinChain_source,
                        proteinChain_target,
                        proteinSequence_source_FromPDBFile,
                        proteinSequence_target_FromPDBFile,
                        nodeName_source,
                        nodeName_target,
                        useAlphaFold);

                    if (string.IsNullOrWhiteSpace(finalStr))
                        Console.WriteLine("ERROR: Error retrieving the PDB file for protein " + ptnSource.proteinID + ".");

                    bw.Write(finalStr);
                }
            }
            catch (IOException ex)
            {
                Console.WriteLine($"ERROR: Problems while writing script.");
                finalStr = "";
            }
            if (string.IsNullOrEmpty(finalStr))
                return "ERROR";

            return filePath.FullName;
        }
        internal static string CreatePythonScript(
            Protein ptnSource,
            Protein ptnTarget,
            List<Crosslink> crossLinks,
            string pdbFile_source,
            string pdbFile_target,
            bool HasMoreThanOneChain_proteinSource,
            bool HasMoreThanOneChain_proteinTarget,
            string proteinChain_source,
            string proteinChain_target,
            string proteinSequence_source_FromPDBFile,
            string proteinSequence_target_FromPDBFile,
            string nodeName_source,
            string nodeName_target,
            bool useAlphaFold)
        {
            //Console.WriteLine("INFO: Getting PDB file name...");

            string[] pdbFilePathName_source = GetPDBFilePathName(pdbFile_source);
            string[] pdbFilePathName_target = !pdbFile_source.Equals(pdbFile_target) ? GetPDBFilePathName(pdbFile_target) : null;
            var sbScript = new StringBuilder();
            //sbScript.Append("cd ").AppendLine(pdbFilePathName_source[0]);
            //sbScript.Append("load ").AppendLine(pdbFilePathName_source[1]);

            //if (!pdbFile_source.Equals(pdbFile_target))
            //    sbScript.Append("load ").AppendLine(pdbFilePathName_target[1]);
            //sbScript.AppendLine("set ignore_case, 0");
            //sbScript.AppendLine($"select chain_{proteinChain_source}, chain {proteinChain_source}");
            //sbScript.AppendLine($"select chain_{proteinChain_target}, chain {proteinChain_target}");
            //sbScript.AppendLine("hide all");

            //Console.WriteLine("INFO: Computing protein offset...");

            int offsetProteinSource = proteinOffsetInPDBSource;
            if (useAlphaFold)
                offsetProteinSource = 0;
            //int offsetProteinSource = ptnSource.fasta.sequence.IndexOf(proteinSequence_source_FromPDBFile);
            //offsetProteinSource = offsetProteinSource > -1 ? offsetProteinSource : 0;
            string selectedResidueItem = "CA";

            //if ((proteinOffsetInPDBSource - 1) == offsetProteinSource)
            //    offsetProteinSource = 0;
            //else if (offsetProteinSource > 0)
            //    offsetProteinSource = (proteinOffsetInPDBSource - 1) - offsetProteinSource;

            int offsetProteinTarget = proteinOffsetInPDBTarget;
            if (useAlphaFold)
                offsetProteinTarget = 0;
            //int offsetProteinTarget = ptnTarget.fasta.sequence.IndexOf(proteinSequence_target_FromPDBFile);
            //offsetProteinTarget = offsetProteinTarget > -1 ? offsetProteinTarget : 0;
            //if ((proteinOffsetInPDBTarget - 1) == offsetProteinTarget)
            //    offsetProteinTarget = 0;
            //else if (offsetProteinTarget > 0)
            //    offsetProteinTarget = (proteinOffsetInPDBTarget - 1) - offsetProteinTarget;

            var sbDistances = new StringBuilder();
            var sbPositions = new StringBuilder();
            int offset_a = -1;
            int offset_b = -1;
            int pos1 = -1;
            int pos2 = -1;

            foreach (var crossLink in crossLinks)
            {
                string chain_a = "";
                string chain_b = "";
                if (crossLink.sourceAccessionNumberProtein.Equals(nodeName_source))
                {
                    chain_a = proteinChain_source;
                    offset_a = offsetProteinSource;
                }
                else
                {
                    chain_a = proteinChain_target;
                    offset_a = offsetProteinTarget;
                }
                if (crossLink.targetAccessionNumberProtein.Equals(nodeName_target))
                {
                    chain_b = proteinChain_target;
                    offset_b = offsetProteinTarget;
                }
                else
                {
                    chain_b = proteinChain_source;
                    offset_b = offsetProteinSource;
                }
                pos1 = crossLink.sourcePosition + offset_a;
                pos2 = crossLink.targetPosition + offset_b;
                sbDistances.AppendLine($"{chain_a}/{pos1}/{selectedResidueItem}, {chain_b}/{pos2}/{selectedResidueItem}");
                sbPositions.AppendLine($"{chain_a}/{pos1}");
                sbPositions.AppendLine($"{chain_b}/{pos2}");
            }

            var distancesList = new HashSet<string>(sbDistances.ToString().Split(new[] { "\r\n" }, StringSplitOptions.RemoveEmptyEntries));
            //var positionList = new HashSet<string>(sbPositions.ToString().Split(new[] { "\r\n" }, StringSplitOptions.RemoveEmptyEntries));
            int countXL = 1;

            //foreach (var dist in distancesList)
            //{
            //    sbScript.AppendLine($"distance xl{countXL}, {dist}");
            //    countXL++;
            //}

            //sbScript.AppendLine($"select a, res {string.Join("+", positionList)};");
            //sbScript.AppendLine("set dash_width, 5");
            //sbScript.AppendLine("set dash_length, 0.1");
            //sbScript.AppendLine("set dash_color, [1.000, 1.000, 0.000]");
            //sbScript.AppendLine("hide label");
            //sbScript.AppendLine("deselect");
            //sbScript.AppendLine("show sticks, a");
            //sbScript.AppendLine($"show cartoon, chain_{proteinChain_source}");
            //sbScript.AppendLine($"color green, chain {proteinChain_source}");
            //sbScript.AppendLine($"show cartoon, chain_{proteinChain_target}");
            //sbScript.AppendLine($"color red, chain {proteinChain_target}");
            //sbScript.AppendLine($"orient chain_{proteinChain_source} chain_{proteinChain_target}");
            //sbScript.AppendLine("deselect");

            sbScript.AppendLine("from pymol import cmd");
            sbScript.AppendLine($"cmd.load(r\"{pdbFile_source}\")");
            if (!pdbFile_source.Equals(pdbFile_target))
                sbScript.AppendLine($"cmd.load(r\"{pdbFile_target}\")");
            sbScript.AppendLine("_list = []");

            //countXL = 1;

            foreach (var dist in distancesList)
            {
                string[] cols_dist = Regex.Split(dist, ",");
                sbScript.AppendLine($"_list.append([\"xl{countXL}\", cmd.distance(\"xl{countXL}\", \"{cols_dist[0]}\", \"{cols_dist[1]}\")])");
                countXL++;
            }
            sbScript.AppendLine("list_string = \"[\" + \", \".join([f\"[{item[0]}, {item[1]}]\" for item in _list]) + \"]\"");
            sbScript.AppendLine("print(list_string)");

            //DirectoryInfo dr = new DirectoryInfo(Path.GetTempPath() + "\\cytoTmpScripts");
            //string prefix = ptnSource.proteinID + "_" + ptnTarget.proteinID;
            //outputfile = Path.Combine(dr.FullName, $"{prefix}_scr_{Guid.NewGuid()}_distances_output.txt");

            //sbScript.AppendLine($"with open('{outputfile}', 'w') as f: f.write(list_string)");

            return sbScript.ToString();
        }

        internal static string CreatePyMOLScript(
            Protein ptnSource,
            Protein ptnTarget,
            List<Crosslink> crossLinks,
            string pdbFile_source,
            string pdbFile_target,
            bool HasMoreThanOneChain_proteinSource,
            bool HasMoreThanOneChain_proteinTarget,
            string proteinChain_source,
            string proteinChain_target,
            string proteinSequence_source_FromPDBFile,
            string proteinSequence_target_FromPDBFile,
            string nodeName_source,
            string nodeName_target,
            bool useAlphaFold)
        {
            //Console.WriteLine("INFO: Getting PDB file name...");

            string[] pdbFilePathName_source = GetPDBFilePathName(pdbFile_source);
            string[] pdbFilePathName_target = !pdbFile_source.Equals(pdbFile_target) ? GetPDBFilePathName(pdbFile_target) : null;
            var sbScript = new StringBuilder();
            sbScript.Append("cd ").AppendLine(pdbFilePathName_source[0]);
            sbScript.Append("load ").AppendLine(pdbFilePathName_source[1]);

            if (!pdbFile_source.Equals(pdbFile_target))
                sbScript.Append("load ").AppendLine(pdbFilePathName_target[1]);
            sbScript.AppendLine("set ignore_case, 0");
            sbScript.AppendLine($"select chain_{proteinChain_source}, chain {proteinChain_source}");
            sbScript.AppendLine($"select chain_{proteinChain_target}, chain {proteinChain_target}");
            sbScript.AppendLine("hide all");

            //Console.WriteLine("INFO: Computing protein offset...");

            int offsetProteinSource = proteinOffsetInPDBSource;
            if (useAlphaFold)
                offsetProteinSource = 0;
            //int offsetProteinSource = ptnSource.fasta.sequence.IndexOf(proteinSequence_source_FromPDBFile);
            //offsetProteinSource = offsetProteinSource > -1 ? offsetProteinSource : 0;
            string selectedResidueItem = "CA";

            //if ((proteinOffsetInPDBSource - 1) == offsetProteinSource)
            //    offsetProteinSource = 0;
            //else if (offsetProteinSource > 0)
            //    offsetProteinSource = (proteinOffsetInPDBSource - 1) - offsetProteinSource;

            int offsetProteinTarget = proteinOffsetInPDBTarget;
            if (useAlphaFold)
                offsetProteinTarget = 0;
            //int offsetProteinTarget = ptnTarget.fasta.sequence.IndexOf(proteinSequence_target_FromPDBFile);
            //offsetProteinTarget = offsetProteinTarget > -1 ? offsetProteinTarget : 0;
            //if ((proteinOffsetInPDBTarget - 1) == offsetProteinTarget)
            //    offsetProteinTarget = 0;
            //else if (offsetProteinTarget > 0)
            //    offsetProteinTarget = (proteinOffsetInPDBTarget - 1) - offsetProteinTarget;

            var sbDistances = new StringBuilder();
            var sbPositions = new StringBuilder();
            int offset_a = -1;
            int offset_b = -1;
            int pos1 = -1;
            int pos2 = -1;

            foreach (var crossLink in crossLinks)
            {
                string chain_a = "";
                string chain_b = "";
                if (crossLink.sourceAccessionNumberProtein.Equals(nodeName_source))
                {
                    chain_a = proteinChain_source;
                    offset_a = offsetProteinSource;
                }
                else
                {
                    chain_a = proteinChain_target;
                    offset_a = offsetProteinTarget;
                }
                if (crossLink.targetAccessionNumberProtein.Equals(nodeName_target))
                {
                    chain_b = proteinChain_target;
                    offset_b = offsetProteinTarget;
                }
                else
                {
                    chain_b = proteinChain_source;
                    offset_b = offsetProteinSource;
                }
                pos1 = crossLink.sourcePosition + offset_a;
                pos2 = crossLink.targetPosition + offset_b;
                sbDistances.AppendLine($"{chain_a}/{pos1}/{selectedResidueItem}, {chain_b}/{pos2}/{selectedResidueItem}");
                sbPositions.AppendLine($"{chain_a}/{pos1}");
                sbPositions.AppendLine($"{chain_b}/{pos2}");
            }

            var distancesList = new HashSet<string>(sbDistances.ToString().Split(new[] { "\r\n" }, StringSplitOptions.RemoveEmptyEntries));
            var positionList = new HashSet<string>(sbPositions.ToString().Split(new[] { "\r\n" }, StringSplitOptions.RemoveEmptyEntries));
            int countXL = 1;

            foreach (var dist in distancesList)
            {
                sbScript.AppendLine($"distance xl{countXL}, {dist}");
                countXL++;
            }

            sbScript.AppendLine($"select a, res {string.Join("+", positionList)};");
            sbScript.AppendLine("set dash_width, 5");
            sbScript.AppendLine("set dash_length, 0.1");
            sbScript.AppendLine("set dash_color, [1.000, 1.000, 0.000]");
            sbScript.AppendLine("hide label");
            sbScript.AppendLine("deselect");
            sbScript.AppendLine("show sticks, a");
            sbScript.AppendLine($"show cartoon, chain_{proteinChain_source}");
            sbScript.AppendLine($"color green, chain {proteinChain_source}");
            sbScript.AppendLine($"show cartoon, chain_{proteinChain_target}");
            sbScript.AppendLine($"color red, chain {proteinChain_target}");
            sbScript.AppendLine($"orient chain_{proteinChain_source} chain_{proteinChain_target}");
            sbScript.AppendLine("deselect");

            return sbScript.ToString();
        }
        private static string[] GetPDBFilePathName(string pdbFile)
        {
            string separator = Path.DirectorySeparatorChar.ToString();
            string[] pdbPathAndFileName = pdbFile.Split(separator);
            StringBuilder pdbFilePath = new StringBuilder();

            for (int i = 0; i < pdbPathAndFileName.Length - 1; i++)
            {
                pdbFilePath.Append(pdbPathAndFileName[i] + separator);
            }

            string pdbFileName = pdbPathAndFileName[pdbPathAndFileName.Length - 1];
            string path = pdbFilePath.ToString().Substring(0, pdbFilePath.Length - 1);

            return new string[] { path, pdbFileName };
        }

        internal static string GetProteinInfoFromUniprot(string proteinID)
        {
            string responseString = "";
            string _url = "https://www.uniprot.org/uniprot/" + proteinID + ".xml";
            HttpWebRequest request = (HttpWebRequest)WebRequest.Create(_url);
            request.Method = "GET";
            request.Timeout = REQUEST_TIMEOUT;
            request.ReadWriteTimeout = REQUEST_TIMEOUT;
            using (HttpWebResponse response = (HttpWebResponse)request.GetResponse())
            {
                if (response.StatusCode == HttpStatusCode.OK)
                {
                    using (Stream responseStream = response.GetResponseStream())
                    {
                        if (responseStream != null)
                        {
                            StreamReader reader = new StreamReader(responseStream);
                            responseString = reader.ReadToEnd();
                        }
                    }
                }
            }
            return responseString;
        }
        internal static FileInfo GetTmpFile(string prefix, string suffix)
        {
            DirectoryInfo dr = new DirectoryInfo(Path.GetTempPath() + "\\cytoTmpScripts");
            if (!dr.Exists)
                dr.Create();
            try
            {
                return new FileInfo(Path.Combine(dr.FullName, $"{prefix}_scr_{Guid.NewGuid()}.{suffix}"));
            }
            catch (IOException e)
            {
                Console.WriteLine("ERROR: Could not work with tmp dir: " + dr.FullName);
            }
            return null;
        }

        internal static string[] GetPDBorCIFfileFromServer(string pdbID)
        {
            try
            {
                string _url = "https://files.rcsb.org/view/" + pdbID + ".pdb";
                HttpWebRequest request = (HttpWebRequest)WebRequest.Create(_url);
                request.Method = "GET";
                request.Timeout = REQUEST_TIMEOUT;
                using (HttpWebResponse response = (HttpWebResponse)request.GetResponse())
                {
                    if (response.StatusCode == HttpStatusCode.OK)
                    {
                        // Get Response
                        Stream inputStream = response.GetResponseStream();
                        if (inputStream != null)
                        {
                            using (StreamReader rd = new StreamReader(inputStream))
                            {
                                StringBuilder responseBody = new StringBuilder();
                                int totalLines = (int)response.ContentLength;
                                int oldProgress = 0;
                                int summaryProcessed = 0;
                                string line;
                                while ((line = rd.ReadLine()) != null)
                                {
                                    responseBody.AppendLine(line);
                                    summaryProcessed += line.Length + 1;
                                    int newProgress = (int)((double)summaryProcessed / totalLines * 100);
                                    if (newProgress > oldProgress)
                                    {
                                        oldProgress = newProgress;
                                        int currentLineCursor = Console.CursorTop;
                                        Console.SetCursorPosition(0, Console.CursorTop);
                                        Console.Write("INFO: Downloading PDB file from server: " + oldProgress + "%");
                                        Console.SetCursorPosition(0, currentLineCursor);
                                    }
                                }
                                return new string[] { "PDB", responseBody.ToString() };
                            }
                        }
                    }
                    else if (response.StatusCode == HttpStatusCode.NotFound)
                    {
                        Console.WriteLine("WARN: There is no PDB for this ID: " + pdbID + ". Trying to retrieve CIF file...");
                        return new string[] { "CIF", GetCiFfileFromServer(pdbID) };
                    }
                    else
                    {
                        Console.WriteLine("WARN: There is no PDB for this ID: " + pdbID + ".");
                        return new string[] { "" };
                    }
                }

                Console.WriteLine("WARN: There is no PDB for this ID: " + pdbID + ".");
                return new string[] { "" };
            }
            catch (Exception e)
            {
                Console.WriteLine("ERROR: " + e.Message);
                return new string[] { "" };
            }
        }

        internal static string GetCiFfileFromServer(string pdbID)
        {
            try
            {
                string _url = "https://files.rcsb.org/view/" + pdbID + ".cif";
                HttpWebRequest request = (HttpWebRequest)WebRequest.Create(_url);
                request.Method = "GET";
                request.Timeout = REQUEST_TIMEOUT;
                using (HttpWebResponse response = (HttpWebResponse)request.GetResponse())
                {
                    if (response.StatusCode == HttpStatusCode.OK)
                    {
                        // Get Response
                        Stream inputStream = response.GetResponseStream();
                        if (inputStream != null)
                        {
                            using (StreamReader rd = new StreamReader(inputStream))
                            {
                                StringBuilder responseBody = new StringBuilder();
                                int totalLines = (int)response.ContentLength;
                                int oldProgress = 0;
                                int summaryProcessed = 0;
                                string line;
                                while ((line = rd.ReadLine()) != null)
                                {
                                    responseBody.AppendLine(line);
                                    summaryProcessed += line.Length + 1;
                                    int newProgress = (int)((double)summaryProcessed / totalLines * 100);
                                    if (newProgress > oldProgress)
                                    {
                                        oldProgress = newProgress;
                                        int currentLineCursor = Console.CursorTop;
                                        Console.SetCursorPosition(0, Console.CursorTop);
                                        Console.Write("INFO: Downloading PDB file from server: " + oldProgress + "%");
                                        Console.SetCursorPosition(0, currentLineCursor);
                                    }
                                }
                                return responseBody.ToString();
                            }
                        }
                    }
                    else
                    {
                        Console.WriteLine("WARN: There is no CIF for this ID: " + pdbID + ".");
                        return "";
                    }
                }
                Console.WriteLine("WARN: There is no CIF for this ID: " + pdbID + ".");
                return "";
            }
            catch (Exception e)
            {
                Console.WriteLine("ERROR: " + e.Message);
                return "";
            }
        }

        internal static string[] GetPDBfileFromSwissModelServer(string _url)
        {
            try
            {
                HttpWebRequest request = (HttpWebRequest)WebRequest.Create(_url);
                request.Method = "GET";
                request.Accept = "text/html";
                request.Headers["Accept-Language"] = "en-US";
                request.Headers["Connection"] = "close";
                request.Timeout = REQUEST_TIMEOUT;
                using (HttpWebResponse response = (HttpWebResponse)request.GetResponse())
                {
                    if (response.StatusCode == HttpStatusCode.OK)
                    {
                        // Get Response
                        Stream inputStream = response.GetResponseStream();
                        if (inputStream != null)
                        {
                            using (StreamReader rd = new StreamReader(inputStream))
                            {
                                StringBuilder responseBody = new StringBuilder();
                                int totalLines = (int)response.ContentLength;
                                int oldProgress = 0;
                                int summaryProcessed = 0;
                                string line;
                                while ((line = rd.ReadLine()) != null)
                                {
                                    responseBody.AppendLine(line);
                                    summaryProcessed += line.Length + 1;
                                    int newProgress = (int)((double)summaryProcessed / totalLines * 100);
                                    if (newProgress > oldProgress)
                                    {
                                        oldProgress = newProgress;
                                        int currentLineCursor = Console.CursorTop;
                                        Console.SetCursorPosition(0, Console.CursorTop);
                                        Console.Write("INFO: Downloading PDB file from server: " + oldProgress + "%");
                                        Console.SetCursorPosition(0, currentLineCursor);
                                    }
                                }
                                return new string[] { "PDB", responseBody.ToString() };
                            }
                        }
                    }
                    else
                    {
                        Console.WriteLine("WARN: There is no PDB file.");
                        return new string[] { "" };
                    }
                }

                Console.WriteLine("WARN: There is no PDB file.");
                return new string[] { "" };
            }
            catch (Exception e)
            {
                Console.WriteLine("ERROR: " + e.Message);
                return new string[] { "" };
            }
        }

        internal static string[] GetPDBfileFromAlphaFoldServer(string pdbID, bool modifyChain)
        {
            try
            {
                string _url = "https://alphafold.ebi.ac.uk/files/AF-" + pdbID + "-F1-model_v2.pdb";
                HttpWebRequest request = (HttpWebRequest)WebRequest.Create(_url);
                request.Method = "GET";
                request.Accept = "text/html";
                request.Headers["Accept-Language"] = "en-US";
                request.Headers["Connection"] = "close";
                request.Timeout = REQUEST_TIMEOUT;
                request.ReadWriteTimeout = REQUEST_TIMEOUT;
                using (HttpWebResponse response = (HttpWebResponse)request.GetResponse())
                using (Stream stream = response.GetResponseStream())
                using (StreamReader reader = new StreamReader(stream))
                {
                    if (response.StatusCode == HttpStatusCode.OK)
                    {
                        string line;
                        StringBuilder responseBuilder = new StringBuilder();
                        int totalLines = (int)response.ContentLength;
                        int oldProgress = 0;
                        int summaryProcessed = 0;
                        while ((line = reader.ReadLine()) != null)
                        {
                            responseBuilder.AppendLine(line);
                            summaryProcessed += line.Length + 1;
                            int newProgress = (int)((double)summaryProcessed / totalLines * 100);
                            if (newProgress > oldProgress)
                            {
                                oldProgress = newProgress;
                                int currentLineCursor = Console.CursorTop;
                                Console.SetCursorPosition(0, Console.CursorTop);
                                Console.Write("INFO: Downloading PDB file from server: " + oldProgress + "%");
                                Console.SetCursorPosition(0, currentLineCursor);
                            }
                        }
                        string finalResponse = responseBuilder.ToString();
                        if (modifyChain)
                        {
                            finalResponse = finalResponse.Replace(" A ", " B ");
                        }
                        return new string[] { "PDB", finalResponse };
                    }
                    else
                    {
                        Console.WriteLine($"INFO: There is no PDB for this ID: {pdbID}.");
                        return new string[] { "" };
                    }
                }
            }
            catch (Exception e)
            {
                Console.WriteLine("ERROR: " + e.Message);
                return new string[] { "" };
            }
        }

        internal static List<FastaItem> GetProteinSequenceFromPDBServer(string pdbID)
        {
            List<FastaItem> fastaList = new();
            try
            {
                string _url = "https://www.rcsb.org/fasta/entry/" + pdbID + "/download";
                HttpWebRequest request = (HttpWebRequest)WebRequest.Create(_url);
                request.Method = "GET";
                request.Accept = "text/html";
                request.Headers["Accept-Language"] = "en-US";
                request.Headers["Connection"] = "close";
                request.Timeout = REQUEST_TIMEOUT;
                using (HttpWebResponse response = (HttpWebResponse)request.GetResponse())
                using (Stream inputStream = response.GetResponseStream())
                using (StreamReader rd = new StreamReader(inputStream))
                {
                    StringBuilder responseData = new StringBuilder();
                    if (response.StatusCode == HttpStatusCode.OK)
                    {
                        string line;
                        while ((line = rd.ReadLine()) != null)
                        {
                            responseData.AppendLine(line);
                        }
                    }
                    else
                    {
                        return new List<FastaItem>();
                    }
                    string[] cols = responseData.ToString().Split('\r');
                    for (int i = 0; i < cols.Length - 1; i += 2)
                    {
                        string seq = cols[i + 1];
                        seq = Regex.Replace(seq, "\n", "").Trim();
                        fastaList.Add(new FastaItem(cols[i].Trim(), seq, seq.Length));
                    }
                }
            }
            catch (Exception e)
            {
                Console.WriteLine("ERROR: Download fasta file from RCSB PDB server:" + e.Message);
                return new List<FastaItem>();
            }
            return fastaList;
        }

        internal static List<FastaItem> GetProteinSequencesFromPDBFile(string fileName, Protein ptn)
        {
            Dictionary<string, int> ResiduesDict = CreateResiduesDict();
            StringBuilder sbSequence = new StringBuilder();
            string proteinChain = "";
            List<FastaItem> fastaList = new();
            try
            {
                using (StreamReader parserFile = new StreamReader(fileName))
                {
                    string line;
                    string lastOffset = "-1";
                    int threshold = 15; // Number of amino acids
                    int countAA = 0;
                    bool getSequence = true;
                    while ((line = parserFile.ReadLine()) != null)
                    {
                        try
                        {
                            if (string.IsNullOrEmpty(line))
                                continue;
                            if (!getSequence)
                            {
                                if (line.StartsWith("ATOM"))
                                    continue;
                                else
                                {
                                    getSequence = true;
                                    countAA = 0;
                                    sbSequence.Clear();
                                }
                            }
                            string[] cols = line.Split(new[] { ' ' }, StringSplitOptions.RemoveEmptyEntries);
                            if (!cols[0].Equals("ATOM"))
                                continue;
                            if (!getSequence)
                            {
                                if (cols[4].Equals(proteinChain))
                                    continue;
                            }
                            proteinChain = cols[4];
                            string pdbResidue = cols[3]; // Residue -> three characters

                            int newResidue = -1;
                            ResiduesDict.TryGetValue(pdbResidue, out newResidue);
                            string currentOffset = cols[5];

                            if (!currentOffset.Equals(lastOffset))
                            {
                                byte[] _byte = new byte[1];
                                _byte[0] = (byte)newResidue;
                                string character = System.Text.Encoding.ASCII.GetString(_byte);
                                sbSequence.Append(character);
                                countAA++;
                                if (countAA > threshold)
                                {
                                    getSequence = false;
                                    // Create Fasta
                                    int proteinOffsetInPDBSource = int.Parse(cols[5]);
                                    string header = ">" + ptn.proteinID + "|Chain " + proteinChain + "|description";
                                    FastaItem fasta = new FastaItem(header, sbSequence.ToString(), sbSequence.ToString().Length, proteinOffsetInPDBSource);
                                    fastaList.Add(fasta);
                                }
                            }
                            lastOffset = currentOffset;
                        }
                        catch (Exception e)
                        {
                            // Handle exception
                        }
                    }
                }
            }
            catch (Exception e)
            {
                Console.WriteLine("ERROR: Problems while reading PDB file: " + fileName);
                return new List<FastaItem>();
            }
            return fastaList;
        }

        internal static Dictionary<string, int> CreateResiduesDict()
        {
            Dictionary<string, int> ResiduesDict = new(); // e.g <GLU, E>
            ResiduesDict.Add("GLY", 71); // Glycine (G)
            ResiduesDict.Add("ALA", 65); // Alanine (A)
            ResiduesDict.Add("SER", 83); // Serine (S)
            ResiduesDict.Add("PRO", 80); // Proline (P)
            ResiduesDict.Add("VAL", 86); // Valine (V)
            ResiduesDict.Add("THR", 84); // Threonine (T)
            ResiduesDict.Add("CYS", 67); // Cystein (C)
            ResiduesDict.Add("ILE", 73); // Isoleucine (I)
            ResiduesDict.Add("LEU", 76); // Leucine (L)
            ResiduesDict.Add("ASN", 78); // Asparagine (N)
            ResiduesDict.Add("ASP", 68); // Aspartic Acid (D)
            ResiduesDict.Add("GLN", 81); // Glutamine (Q)
            ResiduesDict.Add("LYS", 75); // Lysine (K)
            ResiduesDict.Add("GLX", 90); // Glutamic Acid or Glutamine (Z)
            ResiduesDict.Add("GLU", 69); // Glutamic Acid (E)
            ResiduesDict.Add("MET", 77); // Methionine (M)
            ResiduesDict.Add("HIS", 72); // Histidine (H)
            ResiduesDict.Add("PHE", 70); // Phenilanyne (F)
            ResiduesDict.Add("SDC", 85); // Selenocysteine (U)
            ResiduesDict.Add("ARG", 82); // Arginine (R)
            ResiduesDict.Add("TYR", 89); // Tyrosine (Y)
            ResiduesDict.Add("TRP", 87); // Tyrosine (Y)
            return ResiduesDict;
        }

    }
}
