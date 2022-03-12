using System;
using System.Collections.Generic;
using System.Text;

namespace Optimum
{
    public class Graph
    {
        public List<Vertex> allvertexs = new List<Vertex>();
        public List<Edge> alledges = new List<Edge>();

        public void AddVertex(Vertex ver)
        {
            if (allvertexs.Contains(ver) || allvertexs == null) return;
            allvertexs.Add(ver);
        }

        public void AddEdge(Vertex first, Vertex end, double len = 1.0)
        {
            var ed = new Edge(first, end, len);
            if (ed != null/*&& vertexs.Contains(ed.First)&&vertexs.Contains(ed.End)*/)
            {
                alledges.Add(ed);
                if (!first.Edges.Contains(ed) || !end.Edges.Contains(ed))
                {
                    first.Edges.Add(ed);
                    end.Edges.Add(ed);
                }
            }

        }
        public void RemoveVertex(Vertex ver)
        {
            allvertexs.Remove(ver);
        }
        public void Connect(Vertex f, Vertex e, double len = 1.0)
        {
            if (allvertexs.Contains(f) && allvertexs.Contains(e))
            {
                this.AddEdge(f, e);
                return;
            }
            else return;
        }
        public Edge FindEdge(Vertex start, Vertex end)
        {
            foreach(Edge e in alledges)
            {
                if ((start == e.First && end == e.End) || (start == e.End && end == e.First)) return e;
            }
            return null;
        }
        public void ViewGraph()
        {
            for (int i = 0; i < allvertexs.Count; i++)
            {
                allvertexs[i].ViewEdges();
                Console.WriteLine();

            }
            Console.WriteLine();
        }
        public void ViewGraphEdges()
        {
            foreach (var e in alledges)
            {
                Console.WriteLine(e);
            }
        }
        public void BFS(Vertex start)
        {
            //init
            foreach (var v in allvertexs)
            {
                v.Weight = Double.MaxValue;
                v.prev = null;
                v.visited = false;
            }
            start.Weight = 0;
            start.visited = true;
            start.prev = null;
            Queue<Vertex> que = new Queue<Vertex>();
            que.Enqueue(start);//включать De выключать
            while (que.Count > 0)
            {

                Vertex u = que.Dequeue();
                foreach (Edge e in u.Edges)
                {
                    Vertex r = e.End;
                    if (r.visited == false)
                    {
                        r.visited = true;
                        r.Weight = u.Weight + 1;
                        r.prev = u;
                        que.Enqueue(r);
                    }
                }
                u.visited = true;
            }
        }
        public void ViewBFS(Vertex start)
        {
            Console.WriteLine("All Path to {0}", start);
            foreach (var v in allvertexs)
            {
                if (v != start)
                {
                    var tmp = v;
                    while (tmp != null)
                    {
                        Console.Write("{0} => ", tmp); tmp = tmp.prev;
                    }
                    Console.WriteLine();
                }
            }
        }
        public void DFS(Vertex startVertex) //Поиск в глубину (обход графа)
        {
            foreach (Vertex vv in allvertexs)
            {
                vv.visited = false;
                vv.prev = null;
            }

            startVertex.time = 0;

            foreach (Vertex vv in allvertexs)
            {
                if (vv.visited == false)
                {
                    DFS_Visit(startVertex);
                }
            }
        }

        public void DFS_Visit(Vertex u)     //Для DFS
        {
            u.visited = true;
            u.discovered = u.time;
            u.time += 1;    //Счетчик времени увеличиваем на 1

            foreach (Edge ee in u.Edges)
            {
                Vertex v = ee.End;
                if (v.visited == false)
                {
                    v.prev = u;
                    DFS_Visit(v);
                }

            }
        }

        public void PrintPath(Vertex start, Vertex end)
        {
            if (start == end) { Console.WriteLine("{0}", start); }
            else if (end.prev == null) { Console.WriteLine("НЕТ пути "); }
            else PrintPath(start, end.prev); Console.WriteLine(end);

        }
    }
}
