using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace TopologyDistance.Model
{
    public class Crosslink
    {
        public string sourceGene { get; set; }
        public string sourceAccessionNumberProtein { get; set; }
        public string targetGene { get; set; }
        public string targetAccessionNumberProtein { get; set; }
        public int sourcePosition { get; set; }
        public int targetPosition { get; set; }
        public double score { get; set; }
        public double distance { get; set; }

        public Crosslink(string sourceGene, string targetGene, int sourcePosition, int targetPosition, double score)
        {
            this.sourceGene = sourceGene;
            this.targetGene = targetGene;
            this.sourcePosition = sourcePosition;
            this.targetPosition = targetPosition;
            this.score = score;
        }
    }
}
