using System;
using System.Collections.Generic;
using System.Text;

namespace Optimum
{
   public class Vertex
    {
        public string Name { get; set; }
        public double Weight { get; set; }
        public bool visited { get; set; }
        /*public List<Vertex> AllVertexs=new List<Vertex>();
        public List<Vertex> VertexEdget = new List<Vertex>();*/
        //это для графа
        public List<Edge> Edges = new List<Edge>();
        public int Level { get; set; }
        public Vertex prev;

        //DFS
        public int discovered;  //обнаруженная вершина
        public int finished;    //обработанная вершина
        public int time; //метка времени 


        public Vertex(string name, double weight)
        {
            Name = name;
            Weight = weight;
        }
        public Vertex()
        {
            Name = null;
            Weight = Double.NaN;
        }
        public void ConnectWE(Edge e)
        {
            if (e.First == this) Edges.Add(e);
            else throw new Exception("Conn");
        }

        public void ViewEdges()
        {
            Console.WriteLine(this.ToString());
            foreach (var i in Edges)
            {
                Console.Write("{0} , ", i);

            }
            Console.WriteLine();
        }

        public override string ToString()
        {
            return String.Format("({0},{1})", Name, Weight);
        }
    }
}
