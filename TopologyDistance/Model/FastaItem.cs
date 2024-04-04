using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace TopologyDistance.Model
{
    public class FastaItem
    {
        public string accessionNumber { get; set; }
        public string sequence { get; set; }
        public int length { get; set; }
        public int offset { get; set; } = -1;

        public FastaItem(string accessionNumber, string sequence, int length, int offset) : this(accessionNumber, sequence, length)
        {
            this.offset = offset;
        }

        public FastaItem(string accessionNumber, string sequence, int length)
        {
            this.accessionNumber = accessionNumber;
            this.sequence = sequence;
            this.length = length;
        }
    }
}
