using System;
using System.Collections.Generic;
using System.Text;

namespace Optimum
{
   public class Edge
    {
        public Vertex First;
        public Vertex End;
        public double Length { get; set; }
        public Edge(Vertex first, Vertex end, double len = 1.0)
        {
            First = first;
            End = end;
            Length = len;
        }
        public override string ToString()
        {
            return string.Format("Edge ({0}--{1}:{2})", First.ToString(), End.ToString(), Length.ToString());
        }

    }
}
