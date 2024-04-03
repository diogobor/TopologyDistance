using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace TopologyDistance.Model
{
    public class Protein
    {
        public string proteinID { get; set; }
        public string fullName { get; set;}
        public string gene { get; set;}
        public FastaItem fasta { get; set; }
        public string checksum { get; set; }
        public List<PDB> pdbIDs { get; set; }
        public List<Crosslink> monolinks { get; set; }

        public Protein(string proteinID, string fullName, string gene, string sequence, string checksum, List<PDB> pdbIDs, List<Crosslink> monolinks)
        {
            this.proteinID = proteinID;
            this.fullName = fullName;
            this.gene = gene;
            this.fasta = new FastaItem(proteinID, sequence, sequence.Length);
            this.checksum = checksum;
            this.pdbIDs = pdbIDs;
            this.monolinks = monolinks;
        }
        public Protein(string proteinID, string fullName, string gene, FastaItem sequence, string checksum, List<PDB> pdbIDs, List<Crosslink> monolinks)
        {
            this.proteinID = proteinID;
            this.fullName = fullName;
            this.gene = gene;
            this.fasta = sequence;
            this.checksum = checksum;
            this.pdbIDs = pdbIDs;
            this.monolinks = monolinks;
        }

        public Protein(string proteinID, FastaItem sequence, List<Crosslink> monolinks)
        {
            this.proteinID = proteinID;
            this.fasta = sequence;
            this.monolinks = monolinks;
        }

        public Protein()
        {
            pdbIDs = new();
            monolinks = new();
        }
    }
}
