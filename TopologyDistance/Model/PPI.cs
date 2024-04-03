using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace TopologyDistance.Model
{
    public class PPI
    {
        public Protein proteinA { get; set; }
        public Protein proteinB { get; set; }
        public double score { get; set; }
        public List<Crosslink> crosslinks { get; set; }

        public PPI(Protein proteinA, Protein proteinB, double score, List<Crosslink> crosslinks)
        {
            this.proteinA = proteinA;
            this.proteinB = proteinB;
            this.score = score;
            this.crosslinks = crosslinks;
        }
    }
}
