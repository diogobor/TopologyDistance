using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace TopologyDistance.Model
{
    public class PDB
    {
        public string entry { get; set; }
        public string resolution { get; set; }
        public string chain { get; set; }
        public string positions { get; set; }

        public PDB(string entry, string resolution, string chain, string positions)
        {
            this.entry = entry;
            this.resolution = resolution;
            this.chain = chain;
            this.positions = positions;
        }
    }
}
