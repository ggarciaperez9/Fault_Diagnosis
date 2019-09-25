using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;
using System.Linq.Expressions;
using DataAnalysis.Basic;
using DataAnalysis.Statistical;
using GraphVizWrapper;
using GraphVizWrapper.Tests;
using GraphVizWrapper.Commands;
using GraphVizWrapper.Queries;

namespace DiagnosticoDeFallas
{
    internal enum Colors
    {
        White,
        Grey,
        Black
    }
    class Grafo
    {
        public Grafo()
        {
            GrafoCount = 0;
            ListaVertices = new List<Vertice>();
        }
        public List<Vertice> ListaVertices { get; set; }
        public int GrafoCount { get; set; }

        public void AgregarVertice(Vertice a)
        {
            if (GrafoCount < ListaVertices.Count)
            {
                for (int i = 0; i < ListaVertices.Count; i++)
                {
                    if (ListaVertices[i] == null)
                    {
                        ListaVertices[i]=new Vertice(a.Nombre,a.LimInferior,a.LimSuperior,a.Valor);
                        GrafoCount++;
                        break;
                    }
                }
            }
            else
            {
                ListaVertices.Add(new Vertice(a.Nombre, a.LimInferior, a.LimSuperior,a.Valor));
                GrafoCount++;
            }
        }

        public void AgregarArco(string nombre, string nombre1, int signo)
        {
            int pos1 = 0;
            bool entro = false;
            if ((nombre == null) || (nombre1 == null))
            {
                Console.WriteLine("Escribe un nombre...");
                return;
            }
            for (int i = 0; i < ListaVertices.Count; i++)
            {
                if (ListaVertices[i].Nombre == nombre)
                {
                    entro = true;
                    break;
                }
            }
            if (!entro)
            {
                Console.WriteLine("No existe el Vertice " + nombre);
                entro = false;
                return;
            }

            for (int i = 0; i < ListaVertices.Count; i++)
            {
                if (ListaVertices[i].Nombre == nombre1)
                {
                    pos1 = i;
                    entro = true;
                    break;
                }
            }
            if (!entro)
            {
                Console.WriteLine("No existe el Vertice " + nombre1);
                entro = false;
                return;
            }

            for (int i = 0; i < ListaVertices.Count; i++)
            {
                if (ListaVertices[i].Nombre == nombre)
                {
                    ListaVertices[i].AnadirArco(nombre1,pos1,signo);
                    break;
                }
            }
        }

        public void MostrarVertices()
        {
            for (int i = 0; i < ListaVertices.Count; i++)
            {
                if (ListaVertices[i] != null)
                {
                    Console.WriteLine(ListaVertices[i]);
                }
            }
        }

        public void MostrarGrafo()
        {
            for (int i = 0; i < ListaVertices.Count; i++)
            {
                if (ListaVertices[i] != null)
                {
                    Vertice a = ListaVertices[i];
                    Console.WriteLine();
                    Console.WriteLine(a.Nombre + ":");
                    for (int j = 0; j < a.AdyLista.Count; j++)
                    {
                        Console.Write(ListaVertices[a.AdyLista[j].Pos].Nombre + "="+ a.AdyLista[j].Peso + "=> ");
                    }
                }
            }
            Console.WriteLine();
            Console.WriteLine();
        }

        public struct Alcanzados
        {
            public int sgn;
            public string name;

            public Alcanzados(int s, string n)
            {
                sgn = s;
                name = n;
            }

            public override string ToString()
            {
                return (String.Format("{0}, {1}", sgn, name));
            }

            public int Sgn
            {
                get
                {
                    return sgn;
                }
                set
                {
                    sgn = value;
                }
            }

            public string Name
            {
                get
                {
                    return name;
                }
                set
                {
                    name = value;
                }
            }
        }

        public void RecorridoAnchura(Vertice a)
        {
            var alcanzados = new List<Alcanzados>();
            var alc = new Alcanzados();
            int[] d = new int[ListaVertices.Count];
            int[] pi = new int[ListaVertices.Count];
            Colors[] visit = new Colors[ListaVertices.Count];
            for (int i = 0; i < ListaVertices.Count; i++)
            {
                visit[i] = Colors.White;
            }

            int pos = ListaVertices.IndexOf(a);
            visit[pos] = Colors.Grey;
            d[pos] = 0;

            Queue<int> q = new Queue<int>();
            q.Enqueue(pos);

            while (q.Count != 0)
            {
                int tmp = q.Dequeue();
                foreach (var v in ListaVertices[tmp].AdyLista)
                {
                    if (ValorDelSigno(ListaVertices[tmp], ListaVertices[v.Pos], v))
                    {
                        if (visit[v.Pos] == Colors.White)
                        {
                            alc.name= ListaVertices[v.Pos].Nombre;
                            alc.sgn = ListaVertices[v.Pos].signo;
                            alcanzados.Add(alc);
                            visit[v.Pos] = Colors.Grey;
                            d[v.Pos] = d[tmp] + 1;
                            pi[v.Pos] = tmp;
                            q.Enqueue(v.Pos);
                        }
                    }


                }
                visit[tmp] = Colors.Black;
            }
            if (a.Signo != 0)
            {
                Console.WriteLine();
                
                Console.WriteLine("La falla en la variable " + a.Nombre + " con signo : " + a.signo + " se propaga a las variables:");

                for (int i = 0; i < alcanzados.Count; i++)
                {
                    if (alcanzados[i].sgn != 0)
                    {
                        Console.WriteLine("=> " + alcanzados[i].name + " con signo " + alcanzados[i].sgn);
                    }
                }

            }   
        }

        void CrearTabla()
        {
            //Esta funcion le prmite al usuario analizar el comportamiento
            //de la planta ante una falla especifica y como se propaga dicha falla
            Console.Write("Introduzca el Vertice:");
            string name = Console.ReadLine();
            int sgn = 0;
            for (int i = 0; i < ListaVertices.Count; i++)
            {
                if (ListaVertices[i].nombre == name)
                {
                    Console.WriteLine("Introduzca el signo del Vertice:  " + name);
                    sgn = int.Parse(Console.ReadLine());
                    if ((sgn != 0) && ((sgn == 1) || (sgn == -1)))
                    {
                        ListaVertices[i].signo = sgn;
                        RecorridoAnchura(ListaVertices[i]);
                    }
                    else
                        return;
                    
                }
            }
            //Se analiza el aporte de los errores de las otras
            //variables que fueron afectadas
            for (int j = 0; j < ListaVertices.Count; j++)
            {
                if (ListaVertices[j] != null)
                {
                    if ((ListaVertices[j].Signo != 0) && (ListaVertices[j].nombre != name))
                    {
                        RecorridoAnchura(ListaVertices[j]);
                    }
                }

            }
            //Se llevan las variables al estado normal de operacion
            //para posterior analisis
            for (int i = 0; i < ListaVertices.Count; i++)
            {
                ListaVertices[i].Signo = 0;
            }
        }

        public bool ValorDelSigno(Vertice a, Vertice b, Node ab)
        {
            if (a.Signo * ab.Peso> 0)
            {
                if (b.Signo == -1)
                {
                    b.Signo = 0;
                    return false;
                }
                else
                {
                    b.Signo = 1;
                    return true;
                }
            }
            else
            {
                if (b.Signo == 1)
                {
                    b.Signo = 0;
                    return false;
                }
                else
                {
                    b.Signo = -1;
                    return true;
                }
            }
        }

        public void SuprimirVertice(string nombre)
        {
            int pos = 0;
            bool entro = false;
            for (int i = 0; i < ListaVertices.Count; i++)
            {
                if (ListaVertices[i].Nombre == nombre)
                {
                    pos = i;
                    entro = true;
                    break;
                }
            }
            if (!entro)
            {
                Console.WriteLine("No existe el Vertice " + nombre);
                entro = false;
                return;
            }

            for (int i = 0; i < ListaVertices.Count; i++)
            {
                for (int j = 0; j < ListaVertices[i].AdyLista.Count; j++)
                {
                    if (ListaVertices[i].AdyLista[j].Pos == pos)
                    {
                        ListaVertices[i].AdyLista.Remove(ListaVertices[i].AdyLista[j]);
                        break;
                    }
                }
            }
            ListaVertices[pos] = null;
            GrafoCount--;
        }

        public void SuprimirArco(string nombreOrigen, string nombreDestino)
        {
            int posOrigen = 0;
            int posDestino = 0;
            bool entro = false;
            for (int i = 0; i < ListaVertices.Count; i++)
            {
                if (ListaVertices[i].Nombre == nombreOrigen)
                {
                    posOrigen = i;
                    entro = true;
                    break;
                }
            }
            if (!entro)
            {
                Console.WriteLine("No existe el Vertice " + nombreOrigen);
                entro = false;
                return;
            }
                

            for (int i = 0; i < ListaVertices.Count; i++)
            {
                if (ListaVertices[i].Nombre == nombreDestino)
                {
                    posDestino = i;
                    entro = true;
                    break;
                }
            }
            if (!entro)
            {
                Console.WriteLine("No existe el Vertice " + nombreDestino);
                entro = false;
                return;
            }
            for (int i = 0; i < ListaVertices[posOrigen].AdyLista.Count; i++)
            {
                if (ListaVertices[posOrigen].AdyLista[i].Pos == posDestino)
                    ListaVertices[posOrigen].AdyLista.Remove(ListaVertices[posOrigen].AdyLista[i]);
            }
        }

        public void CargarGrafo()
        {
            try
            {
                //using (StreamReader sr = new StreamReader("C:/datos.txt"))
                using (StreamReader sr = new StreamReader(AppDomain.CurrentDomain.BaseDirectory + "/datos_prueba.txt"))
                {
                    //Causality.Visualizar();
                    string lines;
                    var valores = new List<int[]>();
                    lines = sr.ReadLine();
                    var temp = lines.Split('\t');
                    Grafo g = new Grafo();
                    //Se Crean los Nodos
                    for (int i = 0; i < temp.Length; i++)
                    {
                        Vertice a = new Vertice(temp[i]);
                        g.AgregarVertice(a);
                    }
                    var causalMatrix = new int[g.ListaVertices.Count, g.ListaVertices.Count];
                    var rootCausePL = new Dictionary<string, double>();

                    //g.Pearson(causalMatrix, g.GrafoCount);
                    Causality.ObtainCausality();

                    //Se obtiene la Lista de Prioridad de la Causa Raiz
                    double temporal = 0.0;
                    for (int i = 0; i < g.GrafoCount; i++)
                    {
                        for (int j = 0; j < g.GrafoCount; j++)
                        {
                            temporal += causalMatrix[i, j];
                        }
                        rootCausePL.Add(g.ListaVertices[i].Nombre, temporal);
                    }
                    //Se valida el conocimiento del Proceso con los datos
                    for (int i = 0; i < g.ListaVertices.Count; i++)
                    {
                        lines = sr.ReadLine();
                        if (lines == "")
                        {
                            break;
                        }
                        temp = lines.Split('\t');
                        for (int j = 0; j < temp.Length; j++)
                        {
                            for (int k = 0; k < g.ListaVertices.Count; k++)
                            {
                                if (temp[j] == g.ListaVertices[k].Nombre)
                                {
                                    if (causalMatrix[i, k] != 0)
                                    {
                                        g.AgregarArco(g.ListaVertices[i].Nombre, temp[j], causalMatrix[i, k]);
                                    }
                                }
                            }
                        }
                    }
                    //Se ponen los Limites Inferiores
                    lines = sr.ReadLine();
                    temp = lines.Split('\t');
                    for (int j = 0; j < temp.Length; j++)
                    {
                        g.ListaVertices[j].LimInferior = int.Parse(temp[j]);
                    }
                    //Se ponen los Limites Superiores
                    lines = sr.ReadLine();
                    temp = lines.Split('\t');
                    for (int j = 0; j < temp.Length; j++)
                    {
                        g.ListaVertices[j].LimSuperior = int.Parse(temp[j]);
                    }
                    
                    g.MostrarVertices();
                    g.MostrarGrafo();
                    g.CrearTabla();
                    Console.WriteLine();
                    //g.Visualizar();
                }
            }
            catch (Exception e)
            {
                //Console.WriteLine("No se pudo leer el archivo");
                Console.WriteLine(e.Message);
            }
        }

        //private void Visualizar()
        //{
        //    // These three instances can be injected via the IGetStartProcessQuery, 
        //    //                                               IGetProcessStartInfoQuery and 
        //    //                                               IRegisterLayoutPluginCommand interfaces

        //    var getStartProcessQuery = new GetStartProcessQuery();
        //    var getProcessStartInfoQuery = new GetProcessStartInfoQuery();
        //    var registerLayoutPluginCommand = new RegisterLayoutPluginCommand(getProcessStartInfoQuery, getStartProcessQuery);

        //    // GraphGeneration can be injected via the IGraphGeneration interface

        //    var wrapper = new GraphGeneration(getStartProcessQuery,
        //                                      getProcessStartInfoQuery,
        //                                      registerLayoutPluginCommand);

        //    byte[] output = wrapper.GenerateGraph("digraph{a -> b; b -> c; c -> a;}", Enums.GraphReturnType.Pdf);
            
        //}

    }
}

