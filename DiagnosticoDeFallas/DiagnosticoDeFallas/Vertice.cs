using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace DiagnosticoDeFallas
{
    class Vertice
    {
        public string nombre;
        public int limSuperior;
        public int limInferior;
        public double valor;
        public int signo;
        public Vertice(string nombre)
        {
            Nombre = nombre;
            LimSuperior = 0;
            LimInferior = 0;
            Valor = 0;
            Signo = 0;
            AdyLista = new List<Node>();
        }
        public Vertice(string nombre, int limInferior, int limSuperior,double valor)
        {
            Nombre = nombre;
            LimSuperior = limSuperior;
            LimInferior = limInferior;
            Valor = valor;
            if (Valor > LimSuperior)
                Signo = 1;
            else if (Valor < LimInferior)
                Signo = -1;
            else
                Signo = 0;
            AdyLista=new List<Node>();
        }

        public string Nombre
        {
            get { return nombre; }
            set { nombre = value; }
        }
        public int LimSuperior
        {
            get { return limSuperior; }
            set { limSuperior = value; }
        }
        public int LimInferior
        {
            get { return limInferior; }
            set { limInferior = value; }
        }
        public double Valor
        {
            get { return valor; }
            set
            {
                valor = value;
                if (valor > LimSuperior)
                    Signo = 1;
                else if (valor < LimInferior)
                    Signo = -1;
                else
                    Signo = 0;
            }
        }
        public int Signo
        {
            get { return signo; }
            set { signo = value; }
        }
        public List<Node> AdyLista { get; set; }

        public override bool Equals(object obj)
        {
            Vertice aux = obj as Vertice;
            if (aux != null)
                return Nombre == aux.Nombre;
            else
            {
                Console.WriteLine("No se esta comparando con un Vertice");
                return false;
            }

        }
        public override int GetHashCode()
        {
            return base.GetHashCode();
        }

        public void AnadirArco(string nombre,int pos, int signo)
        {
            AdyLista.Add(new Node(nombre,pos,signo));
        }
        public override string ToString()
        {
            return "Name= "+Nombre+" LimInferior= "+LimInferior+ " LimSuperior= "+LimSuperior+" Signo= "+Signo;
        }
    }
}
