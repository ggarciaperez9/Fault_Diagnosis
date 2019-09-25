using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;


namespace DiagnosticoDeFallas
{
    public delegate void MakeGraph();

    class Program
    {
        static void Main(string[] args)
        {
            string opcion;
            Grafo Grafo = new Grafo();
            string nombre;
            string nombre1;
            do
            {
                Menu();
                opcion = Console.ReadLine();
                switch (opcion)
                {
                    case "1":
                        //Grafo = new Grafo();
                        Console.WriteLine("Grafo Creado");
                        break;

                    case "2":
                        Console.Write("Ingrese el Nombre del Vertice: ");
                        nombre = Console.ReadLine();
                        Console.Write("Ingrese el LimInferior: ");
                        int li = TomarEnteroDato();
                        Console.Write("Ingrese el LimSuperior: ");
                        int ls = TomarEnteroDato();
                        Console.Write("Ingrese el valor: ");
                        int valor = TomarEnteroDato();
                        Grafo.AgregarVertice(new Vertice(nombre,li,ls,valor));
                        break;

                    case "3":
                        Console.Write("Ingrese el Nombre del Vertice de Origen: ");
                        nombre = Console.ReadLine();
                        Console.Write("Ingrese el Nombre del Vertice de Destino: ");
                        nombre1 = Console.ReadLine();
                        Console.Write("Ingrese el signo del Arco [1, -1]: ");
                        int signo = TomarEnteroDato();
                        Grafo.AgregarArco(nombre, nombre1, signo);
                        break;

                    case "4":
                        Console.WriteLine("Los Vertices del Grafo son: ");
                        Grafo.MostrarVertices();
                        Console.WriteLine();
                        break;

                    case "5":
                        Console.WriteLine("El Grafo es el siguiente: ");
                        Grafo.MostrarGrafo();
                        Console.WriteLine();
                        break;

                    case "6":
                        Console.Write("Ingrese el Vertice a Eliminar ====> ");
                        nombre = Console.ReadLine();
                        Grafo.SuprimirVertice(nombre);
                        Console.WriteLine();
                        break;

                    case "7":
                        Console.Write("Ingrese el Vertice Origen =====> ");
                        nombre = Console.ReadLine();
                        Console.Write("Ingrese el Vertice Destino =====> ");
                        nombre1 = Console.ReadLine();
                        Grafo.SuprimirArco(nombre, nombre1);
                        Console.WriteLine("\n");
                        break;

                    case "8":
                        Console.WriteLine("Diagnostico Terminado: ");
                        for (int i = 0; i < Grafo.ListaVertices.Count; i++)
                        {
                            if (Grafo.ListaVertices[i] != null)
                            {
                                if (Grafo.ListaVertices[i].Signo != 0)
                                {
                                    Grafo.RecorridoAnchura(Grafo.ListaVertices[i]);
                                }
                            }

                        }
                        break;

                    case "9":
                        Grafo.CargarGrafo();
                        //var inicio = DateTime.Now;
                        //Action a = Causality.ObtainCausality;
                        //var tarea = Task.Factory.StartNew(a);
                        //var fin = DateTime.Now;
                        //Console.WriteLine(fin - inicio);
                        break;

                }

            } while (opcion != "0");

        }
        public static void Menu()
        {
            Console.WriteLine();
            Console.WriteLine("MENU");
            Console.WriteLine("1-Crear Grafo");
            Console.WriteLine("2-Agregar Vertice");
            Console.WriteLine("3-Agregar Arco");
            Console.WriteLine("4-Mostrar Vertices");
            Console.WriteLine("5-Mostrar Grafo");
            Console.WriteLine("6-Eliminar Vertice");
            Console.WriteLine("7-Eliminar Arco");
            Console.WriteLine("8-Diagnosticar");
            Console.WriteLine("9-Cargar Grafo");
            Console.WriteLine("0-Salir");
            Console.WriteLine("Elija una Opcion:");
        }

        public static int TomarEnteroDato()
        {
            string a = Console.ReadLine();
            while (a==null)
            {
                Console.WriteLine("Intoduzca un entero, por favor coopere");
                a = Console.ReadLine();
            }
            return int.Parse(a);
        }
        
        public static bool PerteneceAlGrafo(Grafo grafo, string nombre)
        {
            for (int i = 0; i < grafo.ListaVertices.Count; i++)
            {
                if (grafo.ListaVertices[i].Nombre == nombre)
                    return true;
            }
            return false;
        }
    }
}
