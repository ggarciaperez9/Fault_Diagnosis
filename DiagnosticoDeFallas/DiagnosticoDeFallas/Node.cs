using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace DiagnosticoDeFallas
{
    public class Node
    {
        public Node(string nombre, int pos, int peso)
        {
            Nombre=nombre;
            Peso = peso;
            Pos = pos;
        }

        public string Nombre { get; set; }
        public int Peso { get; set; }
        public int Pos { get; set; }
    }
}
