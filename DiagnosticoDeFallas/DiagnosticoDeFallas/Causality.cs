using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;
using System.IO;
using DataAnalysis.Basic;
using DataAnalysis.Statistical;
using Microsoft.VisualStudio.GraphModel;



namespace DiagnosticoDeFallas
{
    public class Causality
    {
        public static double[,] CM;
        public struct Correlacion
        {
            public Correlacion(int lamb, double coef)
            {
                Lambda = lamb;
                Coeficiente = coef;
            }

            public int Lambda { get; set; }
            public double Coeficiente { get; set; }
        }

        public struct VarParameters
        { 
            public VarParameters(double c, RVector b, RVector e, double s)
            {
                B = new RVector(b);
                C = c;
                Model = new RMatrix();
                Epsilon = new RVector(e);
                Media = c/(1 - B.Sum());
                Sigma = s;
            }

            public VarParameters(double c, RVector b, RVector e, double s, RMatrix m)
            {
                B = new RVector(b);
                C = c;
                Model = new RMatrix(m);
                Epsilon = new RVector(e);
                Media = c / (1 - B.Sum());
                Sigma = s;
            }

            public RMatrix Model { get; set; }
            public double C { get; set; }
            public RVector B { get; set; }
            public RVector Epsilon { get; set; }
            public double Media { get; set; }
            public double Sigma { get; set; }
        }

        public static void Visualizar(RMatrix Q, List<string> variables, string metodo)
        {
            Graph g = new Graph();
            GraphPropertyCollection properties = g.DocumentSchema.Properties;
            //GraphProperty background = properties.AddNewProperty("Background", typeof(Brush));
            GraphProperty size = properties.AddNewProperty("Size", typeof(String));
            //GraphProperty start = properties.AddNewProperty("Start", typeof(DateTime));

            //GraphNode nodeA = g.Nodes.GetOrCreate("a");
            //nodeA.Label = "Gabriel";
            //nodeA[size] = "10";

            //GraphNode nodeB = g.Nodes.GetOrCreate("b");
            //nodeB.Label = "Danyer";
            ////nodeB[background] = new SolidColorBrush(Color.FromArgb(0xFF, 0x00, 0x80, 0x80));

            //GraphNode nodeC = g.Nodes.GetOrCreate("c");
            //nodeC.Label = "El Compa";
            ////nodeC[start] = new DateTime(2010, 6, 10);

            //// links            
            //g.Links.GetOrCreate(nodeA, nodeB);
            //g.Links.GetOrCreate(nodeA, nodeC);
            var n = new List<GraphNode>();
            for (var i = 0; i < variables.Count; i++)
            {
                n.Add(g.Nodes.GetOrCreate(variables[i]));
                n[i].Label = variables[i];
                n[i][size] = "10";
                
            }
            for (var i = 0; i < Q.NCols; i++)
            {
                for (var j = 0; j < Q.NCols; j++)
                {
                    if (i == j) continue;//Si es la misma variable continua
                    if (Equals(Q[i, j], 0.0)) continue;//Si no existe causalidad continua
                    //if (Q[i, j] < Q[j, i]) continue;//Si existe un lazo se queda el arco mas fuerte
                    g.Links.GetOrCreate(n[i], n[j], Math.Round(Q[i, j], 3).ToString(), null);
                }
            }

            g.Save(AppDomain.CurrentDomain.BaseDirectory + "/" + metodo + ".dgml");
        }

        public static RMatrix TransformRow(RVector v1)
        {
            var result = new RMatrix(1, v1.Length);
            for (int i = 0; i < v1.Length; i++)
                result[0, i] = v1[i];
            return result;
        }

        public static RMatrix TransformCol(RVector v1)
        {
            var result = new RMatrix(v1.Length, 1);
            for (int i = 0; i < v1.Length; i++)
                result[i, 0] = v1[i];
            return result;
        }

        public static RVector Transform(RVector vec, RMatrix mat)
        {
            var result = new RVector(vec.Length);
           if (mat.NRows != vec.Length)
                throw new ArgumentException("The ndim of the vector must be equal"
                                    + " to the number of rows of the matrix!");

            for (int i = 0; i < mat.NRows; i++)
            {
                result[i] = mat[i, 0];
            }
            return result;
        }

        public static RMatrix RemoveRow(int rowToRemove, RMatrix originalArray)
        {
            if ((rowToRemove > originalArray.NRows - 1) || (rowToRemove < 0))
            {
                throw new Exception("El indice de la fila a remover esta fuera de rango");
            }
            var result = new double[originalArray.NRows - 1, originalArray.NCols];

            for (int i = 0, j = 0; i < originalArray.NCols; i++)
            {
                for (int k = 0, u = 0; k < originalArray.NRows; k++)
                {
                    if (k == rowToRemove)
                        continue;

                    result[u, j] = originalArray[k, i];
                    u++;
                }
                j++;
            }

            return new RMatrix(result);
        }

        public static RMatrix RemoveColumn(int columnToRemove, RMatrix originalArray)
        {
            if ((columnToRemove > originalArray.NCols - 1) || (columnToRemove < 0))
            {
                throw new Exception("El indice de la columna a remover esta fuera de rango");
            }
            var result = new double[originalArray.NRows, originalArray.NCols - 1];

            for (int i = 0, j = 0; i < originalArray.NRows; i++)
            {
                for (int k = 0, u = 0; k < originalArray.NCols; k++)
                {
                    if (k == columnToRemove)
                        continue;

                    result[j, u] = originalArray[i, k];
                    u++;
                }
                j++;
            }

            return new RMatrix(result);
        }

        public static bool LeastSquaresTest(RVector y, RVector x)
        {
            var b1 = 0.0;
            var b0 = 0.0;
            var RSS = 0.0;
            var b1Var = 0.0;
            var b0Var = 0.0;
            var b1Error = 0.0;
            var b0Error = 0.0;
            var Sigma = 0.0;
            //var AIC = 0;
            b1 = Summary.Covariance(y, x) / Summary.Variance(x);
            b0 = Summary.Mean(y) - b1 * Summary.Mean(x);

            for (var i = 0; i < x.Length; i++)
                RSS += Math.Pow(y[i] - b0 - b1 * x[i], 2);

            Sigma = RSS/(x.Length - 2); 
            b0Var = Sigma*(1.0/x.Length + Math.Pow(Summary.Mean(x), 2)/(Summary.Variance(x)*x.Length));
            b1Var = Sigma/(Summary.Variance(x)*x.Length);

            b0Error = Math.Sqrt(b0Var);
            b1Error = Math.Sqrt(b1Var);

            var t0 = b0/b0Error;
            var t1 = b1/b1Error;
            //AIC = Causality.AIC(RSS, 2, x.Length);

            return Math.Abs(t1) >= 2.042 ? true : false;
        }

        public static RMatrix OrdinalLeastSquares(RMatrix Y, RMatrix X)
        {
            var bAux = new RMatrix(X.NCols + 1, 1);
            var U = new RVector(X.NRows);
            var RSS = new RVector(X.NRows);
            for (int i = 0; i < U.Length; i++)
            {
                U[i] = 1;
            }
            X.InsertColVector(U, 0);
            bAux = RMatrix.Inverse(X.GetTranspose() * X) * X.GetTranspose() * Y;

            RSS = Transform(RSS, Y - X * bAux);
            var FM = new double();
            for (var i = 0; i < RSS.Length; i++)
            {
                FM += Math.Pow(RSS[i], 2);
            }
            var RM = new double();
            for (var i = 0; i < Y.NRows; i++)
            {
                RM += Math.Pow(Y[i, 0] - Summary.Mean(Y.GetColVector(0)), 2);
            }
            var F = (RM - FM)/(bAux.NRows - 1)/(FM/(Y.NRows - bAux.NRows));

            //var result = new Dictionary<RMatrix, bool>();
            //if (!result.ContainsKey(bAux)) result.Add(bAux, F);

            return bAux;
        }

        public static int UnivarModelOrder(IReadOnlyList<double> y)
        {
            if (y == null) throw new ArgumentNullException(nameof(y));
            var aic = new double[5];
            var RSS = new double[5];
            //var rss = 0.0;
            var N = 0;

            for (var n = 1; n < 6; n++)
            {
                var yActual = new List<double>(y);
                var yPasada = new RMatrix(y.Count - n, n);
                var Y = new RMatrix();
                var k = 0;

                for (var i = n - 1; i >= 0; i--)
                {
                    for (var j = 0; j < y.Count - n; j++)
                        yPasada[j, i] = y[j + k];
                    k++;
                    yActual.RemoveAt(0);
                }
                Y.InsertColVector(new RVector(yActual.ToArray()), 0);

                var I = new RVector(yPasada.NRows);
                var temporal = new RVector(yPasada.NRows);
                for (int i = 0; i < I.Length; i++)
                {
                    I[i] = 1;
                }
                yPasada.InsertColVector(I, 0);
                var bAux = RMatrix.Inverse(yPasada.GetTranspose() * yPasada) * yPasada.GetTranspose() * Y;

                temporal = Transform(temporal, Y - yPasada * bAux);
                for (var m = 0; m < temporal.Length; m++)
                {
                    RSS[n - 1] += Math.Pow(temporal[m], 2);
                }
                aic[n - 1] = Math.Log(RSS[n - 1]) + (double)2 * n / Y.NRows;
            }
            for (var i = 1; i < aic.Length; i++)
            {
                if (aic[i] < aic[N])
                {
                    N = i;
                }
            }
            N = N + 1;
            
            return N;
        }

        public static int MultiModelOrder(int yIndex, RMatrix X)
        {
            var total = 6;
            var aic = new double[total];
            var RSS = new double[total];
            var rss = 0.0;
            var N = 0;

            for (var n = 1; n < total + 1; n++)
            {
                var xAux = new RMatrix();
                var yAux = X.GetColVector(yIndex);
                for (int i = 0, j = 0; i < X.NCols; i++)
                {
                    if (i != yIndex)
                    {
                        xAux.InsertColVector(X.GetColVector(i), j);
                        j++;
                    }

                }
                xAux.InsertColVector(yAux, 0);
                var xActual = new List<double>(xAux.GetColVector(0).ToList());
                var unrestrictedModel = new RMatrix();

                for (var i = 0; i < xAux.NCols; i++)
                {
                    var xPasada = new RMatrix(xAux.GetColVector(i).Length - n, n);
                    var k = 0;

                    for (var m = n - 1; m >= 0; m--)
                    {
                        for (var j = 0; j < xAux.GetColVector(i).Length - n; j++)
                            xPasada[j, m] = xAux.GetColVector(i)[j + k];
                        k++;
                        if (i == 0)
                        {
                            xActual.RemoveAt(0);
                        }
                    }
                    for (var j = 0; j < xPasada.NCols; j++)
                    {
                        unrestrictedModel.InsertColVector(xPasada.GetColVector(j), unrestrictedModel.NCols);
                        
                    }
                }
                var I = new RVector(unrestrictedModel.NRows);

                for (var i = 0; i < I.Length; i++)
                {
                    I[i] = 1;
                }
                var Y = new RMatrix();
                Y.InsertColVector(new RVector(xActual.ToArray()), 0);
                unrestrictedModel.InsertColVector(I, 0);

                var bAux = RMatrix.Inverse(unrestrictedModel.GetTranspose() * unrestrictedModel) * unrestrictedModel.GetTranspose() * Y;

                var temporal = new RVector(unrestrictedModel.NRows);
                temporal = Transform(temporal, Y - unrestrictedModel * bAux);
                
                for (var i = 0; i < temporal.Length; i++)
                {
                    RSS[n - 1] += Math.Pow(temporal[i], 2);
                }
                aic[n - 1] = Math.Log(RSS[n - 1]) + (double)2 * n / Y.NRows;
            }
            for (var i = 1; i < aic.Length; i++)
            {
                if (aic[i] < aic[N])
                {
                    N = i;
                }
            }
            N = N + 1;

            return N;
        }

        public static int[] ModelOrderGC(int yIndex, int uIndex, RMatrix X)
        {
            var rssResModel = new double[6];
            var rssUnResModel = new double[6];
            var aicResModel = new double[6];
            var aicUnResModel = new double[6];
            var nResModel = 0;
            var nUnResModel = 0;

            

            for (var n = 1; n < 6 + 1; n++)
            {
                //var B = new RMatrix();
                var yAux = X.GetColVector(yIndex);
                X = RemoveColumn(yIndex, X);
                X.InsertColVector(yAux, 0);
                var restrictedModel = new RMatrix();
                var unrestrictedModel = new RMatrix();

                for (var i = 0; i < X.NCols; i++)
                {
                    var xActual = new List<double>(X.GetColVector(i).ToList());
                    var xPasada = new RMatrix(X.GetColVector(i).Length - n, n);
                    var k = 0;

                    for (var m = n - 1; m >= 0; m--)
                    {
                        for (var j = 0; j < X.GetColVector(i).Length - n; j++)
                            xPasada[j, m] = X.GetColVector(i)[j + k];
                        k++;
                        xActual.RemoveAt(0);
                    }

                    xPasada.InsertColVector(new RVector(xActual.ToArray()), xPasada.NCols);
                    for (var j = 0; j < xPasada.NCols; j++)
                    {
                        unrestrictedModel.InsertColVector(xPasada.GetColVector(j), unrestrictedModel.NCols);
                        if (i != uIndex)
                        {
                            restrictedModel.InsertColVector(xPasada.GetColVector(j), restrictedModel.NCols);
                        }
                    }
                }
                var I = new RVector(restrictedModel.NRows);
                var temporal = new RVector(restrictedModel.NRows);
                for (var i = 0; i < I.Length; i++)
                {
                    I[i] = 1;
                }
                var Yr = new RMatrix();
                Yr.InsertColVector(restrictedModel.GetColVector(0), 0);
                restrictedModel = RemoveColumn(0, restrictedModel);
                restrictedModel.InsertColVector(I, 0);

                var bAuxR = RMatrix.Inverse(restrictedModel.GetTranspose() * restrictedModel) * restrictedModel.GetTranspose() * Yr;

                temporal = Transform(temporal, Yr - restrictedModel * bAuxR);
                
                for (var i = 0; i < temporal.Length; i++)
                {
                    rssResModel[n - 1] += Math.Pow(temporal[i], 2);
                }
                aicResModel[n - 1] = Math.Log(rssResModel[n - 1]) + 2 * n * Math.Pow(X.NCols - 2, 2)/ Yr.NRows;

                var Yu = new RMatrix();
                Yu.InsertColVector(unrestrictedModel.GetColVector(0), 0);
                unrestrictedModel = RemoveColumn(0, unrestrictedModel);
                unrestrictedModel.InsertColVector(I, 0);

                var bAuxU = RMatrix.Inverse(unrestrictedModel.GetTranspose() * unrestrictedModel) * unrestrictedModel.GetTranspose() * Yu;

                temporal = Transform(temporal, Yu - unrestrictedModel * bAuxU);

                for (var i = 0; i < temporal.Length; i++)
                {
                    rssUnResModel[n - 1] += Math.Pow(temporal[i], 2);
                }
                aicUnResModel[n - 1] = Math.Log(rssUnResModel[n - 1]) + 2 * n * Math.Pow(X.NCols - 1, 2) / Yu.NRows;
            }
            for (var i = 1; i < aicResModel.Length; i++)
            {
                if (aicResModel[i] < aicResModel[nResModel])
                {
                    nResModel = i;
                }
            }
            nResModel = nResModel + 1;

            for (var i = 1; i < aicUnResModel.Length; i++)
            {
                if (aicUnResModel[i] < aicUnResModel[nUnResModel])
                {
                    nUnResModel = i;
                }
            }
            nUnResModel = nUnResModel + 1;

            return new int[]{nResModel, nUnResModel};
        }

        public static int[] Aikake(double[] RM, double[] UM, int K, int p)
        {
            var aicUnResModel = new double[UM.Length];
            var aicResModel = new double[RM.Length];
            var nUnResModel = 0;
            var nResModel = 0;
            for (int i = 0; i < aicUnResModel.Length; i++)
            {
                aicUnResModel[i] = Math.Log(UM[i]) + 2 * (i + 1) * Math.Pow(p - 1, 2) / K;
                aicResModel[i] = Math.Log(RM[i]) + 2 * (i + 1) * Math.Pow(p - 2, 2) / K;
            }

            for (var i = 1; i < aicUnResModel.Length; i++)
            {
                if (aicUnResModel[i] < aicUnResModel[nUnResModel])
                {
                    nUnResModel = i;
                }
            }
            nUnResModel = nUnResModel + 1;

            for (var i = 1; i < aicResModel.Length; i++)
            {
                if (aicResModel[i] < aicResModel[nResModel])
                {
                    nResModel = i;
                }
            }
            nResModel = nResModel + 1;

            return new int[] {nResModel, nUnResModel};
        }

        public static bool F_Test(double ERM, double EURM, int modelOrder, int K, int p)
        {
            var F = ((ERM - EURM)*(K - modelOrder - p*modelOrder)/(modelOrder*EURM));

            return F >= 0.01/(p * (p - 1)) ? true : false;
        }

        //public static RMatrix GrangerCausality(RMatrix X)
        //{
        //    var GC = new double[X.NCols - 1, 7];
        //    var causalMatrix = new RMatrix(X.NCols, X.NCols);
        //    for (var i = 0; i < X.NCols; i++)
        //    {
        //        for (int j = 0, u = 0; j < X.NCols; j++)
        //        {
        //            if (i != j)
        //            {
        //                //var ERM = VarResModel(i, j, X, n[0]);
        //                //var EURM = VarUnResModel(i, j, X, n[1]);
        //                //causalMatrix[i, j] = (1 - EURM/ERM);
        //                //var n = ModelOrderGC(i, j, X);

        //                //var ERM = new double[7];
        //                //var EURM = new double[7];
        //                //for (var k = 1; k < 8; k++)
        //                //{
        //                //    ERM[k - 1] = VarResModel(i, j, X, k);
        //                //    EURM[k - 1] = VarUnResModel(j, X, k);
        //                //}
        //                //var n = Aikake(ERM, EURM, X.GetColVector(j).Length, X.NCols);
        //                //if ((F_Test(ERM[n[0] - 1], EURM[n[1] - 1], n[0], X.GetColVector(j).Length, X.NCols - 1))
        //                //    && (F_Test(ERM[n[0] - 1], EURM[n[1] - 1], n[1], X.GetColVector(j).Length, X.NCols)))
        //                //    causalMatrix[i, j] = (1 - EURM[n[1] - 1] / ERM[n[0] - 1]);
        //                //else
        //                //    causalMatrix[i, j] = 0;

        //                for (var k = 1; k < 8; k++)
        //                {
        //                    //var ERM = VarResModel(i, j, X, k);
        //                    //var EURM = VarUnResModel(j, X, k);
        //                    //GC[u, k - 1] = (1 - EURM / ERM);
        //                    //if (F_Test(ERM, EURM, k, X.GetColVector(j).Length, X.NCols))
        //                    //    GC[u, k - 1] = (1 - EURM / ERM);
        //                    //else
        //                    //    GC[u, k - 1] = 0;

        //                    var ERM = VarResModel(i, j, X, k);
        //                    var EURM = VarUnResModel(j, X, k);
        //                    GC[u, k - 1] = Math.Log(ERM/EURM);
        //                }
        //                u++;
        //            }
        //        }
        //    }

        //    return causalMatrix;
        //}

        //public static RMatrix GCausalityxy(RMatrix X)
        //{
        //    var causalMatrix = new RMatrix(X.NCols, X.NCols);
        //    //var DependVariable = new RMatrix(X.NCols, X.NCols);
        //    for (int i = 0; i < X.NCols; i++)
        //    {
        //        for (int j = 0; j < X.NCols; j++)
        //        {
        //            if (i != j)
        //            {
        //                var M = new RMatrix();
        //                M.InsertColVector(X.GetColVector(j), 0);
        //                M.InsertColVector(X.GetColVector(i), 1);
        //                //M.InsertColVector(X.GetColVector(k), 2);

        //                var nU = MultiModelOrder(0, M);
        //                var EURM = VarUnResModel(0, M, nU);

        //                var Mr = RemoveColumn(1, M);
        //                var nR = MultiModelOrder(0, Mr);
        //                var ERM = VarResModel(1, 0, M, nR);

        //                if (Equals(ERM, 0.0) || Equals(EURM, 0.0))
        //                {
        //                    causalMatrix[i, j] = 0;
        //                    continue;
        //                }
        //                var temp = F_Test(ERM, EURM, nU, X.GetColVector(j).Length, M.NCols) ? 1 - EURM / ERM : 0;

        //                if (causalMatrix[i, j] < temp)
        //                {
        //                    causalMatrix[i, j] = temp;
        //                    //DependVariable[i, j] = k;
        //                }
        //            }
        //        }
        //    }

        //    return causalMatrix;
        //}

        public static RVector GCausality(RMatrix X)
        {
            var causalVector = new RVector(X.NCols);
            var orden = 6;
            var ETA = new RVector(orden);
            var BIC = new RVector(orden);
            for (int n = 0; n < 6; n++)
            {
                var eta = new RMatrix(X.NCols, X.NCols);
                for (int i = 0; i < X.NCols; i++)
                {
                    for (int j = 0; j < X.NCols; j++)
                    {
                        var EURM = VarUnResModel(j, X, n + 1);
                        if (i != j)
                        {
                            //var M = new RMatrix(X);
                            //var nu = MultiModelOrder(j, M);


                            //M = RemoveColumn(i, M);
                            //var nr = MultiModelOrder(j < i ? j : j - 1, M);
                            //var ERM = VarResModel(i, j, X, n);
                            var EURMc = VarUnResModel(i, X, n + 1);
                            var covariance = Summary.Covariance(EURMc, EURM);
                            //if (Equals(ERM, 0.0) || Equals(EURM, 0.0))
                            //{
                            //    causalMatrix[i, j] = 0;
                            //    continue;
                            //}
                            //var temp = F_Test(ERM, EURM, n, X.GetColVector(j).Length, X.NCols) ? Math.Log(ERM/EURM) : 0;
                            eta[i, j] = covariance;
                            //causalMatrix[i, j] = temp;
                            //n[nU - 1, nR - 1] = temp;
                        }
                        else
                        {
                            eta[i, j] = Summary.Variance(EURM);
                        }
                    }
                }
                ETA[n] = RMatrix.LUDeterminant(eta);
            }

            for (int i = 0; i < orden; i++)
            {
                BIC[i] = Math.Log(ETA[i]) + 2*i*Math.Pow(X.NCols, 2)/X.NRows;
            }
            var min = BIC.MinValue();
            for (int i = 0; i < BIC.Length; i++)
            {
                if (Equals(BIC[i], min))
                {
                    for (int j = 0; j < causalVector.Length; j++)
                    {
                        if (j != 0)
                        {
                            var U = Summary.Variance(VarUnResModel(0, X, i + 1));
                            var R = Summary.Variance(VarResModel(j, 0, X, i + 1));
                            causalVector[j] = F_Test(R, U, i + 1, X.GetColVector(j).Length
                                , X.NCols) ? Math.Log(R / U) : 0;
                        }
                        else
                        {
                            causalVector[j] = 0.0;
                        }
                    }
                }
            }
            return causalVector;
        }

        public static RMatrix GCNew(RMatrix X)
        {
            var CausalMatrix = new RMatrix();
            for (int i = 0; i < X.NCols; i++)
            {
                var M = X;
                var t = M.GetColVector(i);
                M = RemoveColumn(i, M);
                M.InsertColVector(t, 0);
                CausalMatrix.InsertRowVector(GCausality(M), i);
            }

            for (int j = 1; j < X.NCols; j++)
            {
                var temp = CausalMatrix.GetRowVector(j - 1);
                for (int k = 0; k < j; k++)
                {
                    var t = temp[k];
                    temp[k - 1] = t;
                }
                temp[j] = 0.0;
            }

            return CausalMatrix;
        }

        public static RVector VarResModel(int uIndex, int yIndex, RMatrix X, int modelOrder)
        {
            //var B = new RMatrix();
            //var modelOrder = ModelOrderGC(yIndex, uIndex, xAux)[0];
            var xAux = new RMatrix();
            var yAux = X.GetColVector(yIndex);
            for (int i = 0, j = 0; i < X.NCols; i++)
            {
                if ((i != yIndex) && (i != uIndex))
                {
                    if (LeastSquaresTest(yAux, X.GetColVector(i)))
                    {
                        xAux.InsertColVector(X.GetColVector(i), j);
                        j++;
                    }
                }
                
            }
            
            //xAux = RemoveColumn(uIndex, xAux);
            //yIndex--;
            
            //xAux = RemoveColumn(yIndex, xAux);
            xAux.InsertColVector(yAux, 0);
            var xActual = new List<double>(xAux.GetColVector(0).ToList());
            var restrictedModel = new RMatrix();
            //var unrestrictedModel = new RMatrix();
            //modelOrder++;
            
            for (var i = 0; i < xAux.NCols; i++)
            {
                //var xActual = new List<double>(xAux.GetColVector(i).ToList());
                var xPasada = new RMatrix(xAux.GetColVector(i).Length - modelOrder, modelOrder);
                var k = 0;

                for (var m = modelOrder - 1; m >= 0; m--)
                {
                    for (var j = 0; j < xAux.GetColVector(i).Length - modelOrder; j++)
                        xPasada[j, m] = xAux.GetColVector(i)[j + k];
                    k++;
                    //xActual.RemoveAt(0);
                    if (i == 0)
                    {
                        xActual.RemoveAt(0);
                    }
                }
                //xPasada.InsertColVector(new RVector(xActual.ToArray()), xPasada.NCols);

                for (var j = 0; j < xPasada.NCols; j++)
                {
                    restrictedModel.InsertColVector(xPasada.GetColVector(j), restrictedModel.NCols);
                }
            }
            var yPasada = new RMatrix();
            for (var i = 0; i < modelOrder; i++)
            {
                yPasada.InsertColVector(restrictedModel.GetColVector(0), i);
                restrictedModel = RemoveColumn(0, restrictedModel);
            }
            var I = new RVector(restrictedModel.NRows);
            
            for (var i = 0; i < I.Length; i++)
                I[i] = 1;

            var Y = new RMatrix();
            Y.InsertColVector(new RVector(xActual.ToArray()), 0);
            restrictedModel.InsertColVector(I, 0);
            
            var bAux = RMatrix.Inverse(restrictedModel.GetTranspose() * restrictedModel) * restrictedModel.GetTranspose() * Y;
            var aAux = RMatrix.Inverse(yPasada.GetTranspose() * yPasada) * yPasada.GetTranspose() * Y;

            var temporal = new RVector(restrictedModel.NRows);
            temporal = Transform(temporal, Y - (yPasada * aAux) - (restrictedModel * bAux));

            var rssResModel = 0.0;
            for (var i = 0; i < temporal.Length; i++)
            {
                rssResModel += Math.Pow(temporal[i], 2);
            }
            //return DurbinWatsonTest(temporal) ? Summary.Variance(temporal) : 0.0;
            //return DurbinWatsonTest(temporal) ? rssResModel : 0.0;
            //var RM = new double();
            //for (var i = 0; i < Y.NRows; i++)
            //{
            //    RM += Math.Pow(Y[i, 0] - Summary.Mean(Y.GetColVector(0)), 2);
            //}
            //var F = (RM - rssResModel) / (bAux.NRows - 1) / (rssResModel / (Y.NRows - bAux.NRows));

            //return rssResModel;
            return temporal;
        }

        public static bool DurbinWatsonTest(RVector e)
        {
            var d = 0.0;
            var sum1 = 0.0;
            var sum2 = 0.0;
            for (int i = 1; i < e.Length; i++)
            {
                sum1 += Math.Pow(e[i] - e[i - 1], 2);
                sum2 += Math.Pow(e[i], 2);
            }
            d = sum1/sum2;
            return 2 <= d && d <= 3;
        }

        public static RVector VarUnResModel(int yIndex, RMatrix X, int modelOrder)
        {
            var xAux = new RMatrix();
            var yAux = X.GetColVector(yIndex);
            for (int i = 0, j = 0; i < X.NCols; i++)
            {
                if (i != yIndex)
                {
                    //if (LeastSquaresTest(yAux, X.GetColVector(i)))
                    //{

                    //}
                    xAux.InsertColVector(X.GetColVector(i), j);
                    j++;
                }

            }
            //var modelOrder = ModelOrderGC(yIndex, uIndex, xAux)[1];
            
            //xAux = RemoveColumn(yIndex, xAux);
            xAux.InsertColVector(yAux, 0);
            var xActual = new List<double>(xAux.GetColVector(0).ToList());
            //var restrictedModel = new RMatrix();
            var unrestrictedModel = new RMatrix();

            for (var i = 0; i < xAux.NCols; i++)
            {
                //var xActual = new List<double>(xAux.GetColVector(i).ToList());
                var xPasada = new RMatrix(xAux.GetColVector(i).Length - modelOrder, modelOrder);
                var k = 0;

                for (var m = modelOrder - 1; m >= 0; m--)
                {
                    for (var j = 0; j < xAux.GetColVector(i).Length - modelOrder; j++)
                        xPasada[j, m] = xAux.GetColVector(i)[j + k];
                    k++;
                    //xActual.RemoveAt(0);
                    if (i == 0)
                    {
                        xActual.RemoveAt(0);
                    }
                }
                //Xt.InsertColVector(new RVector(xActual.ToArray()), 0);
                //xPasada.InsertColVector(new RVector(xActual.ToArray()), xPasada.NCols);
                for (var j = 0; j < xPasada.NCols; j++)
                {
                    unrestrictedModel.InsertColVector(xPasada.GetColVector(j), unrestrictedModel.NCols);
                    //if (i != uIndex)
                    //{
                    //    restrictedModel.InsertColVector(xPasada.GetColVector(j), restrictedModel.NCols);
                    //}
                }
            }
            var yPasada = new RMatrix();
            for (var i = 0; i < modelOrder; i++)
            {
                yPasada.InsertColVector(unrestrictedModel.GetColVector(0), i);
                unrestrictedModel = RemoveColumn(0, unrestrictedModel);
            }
            var I = new RVector(unrestrictedModel.NRows);
            
            for (var i = 0; i < I.Length; i++)
            {
                I[i] = 1;
            }
            var Y = new RMatrix();
            Y.InsertColVector(new RVector(xActual.ToArray()), 0);
            yPasada.InsertColVector(I, 0);

            var bAux = RMatrix.Inverse(unrestrictedModel.GetTranspose() * unrestrictedModel) * unrestrictedModel.GetTranspose() * Y;
            var aAux = RMatrix.Inverse(yPasada.GetTranspose() * yPasada) * yPasada.GetTranspose() * Y;

            var temporal = new RVector(unrestrictedModel.NRows);
            temporal = Transform(temporal, Y - (yPasada * aAux) - (unrestrictedModel * bAux));

            var rssUnResModel = 0.0;
            for (var i = 0; i < temporal.Length; i++)
            {
                rssUnResModel += Math.Pow(temporal[i], 2);
            }

            //return DurbinWatsonTest(temporal) ? Summary.Variance(temporal) : 0.0;
            //return DurbinWatsonTest(temporal) ? rssUnResModel : 0.0;
            //var RM = new double();
            //for (var i = 0; i < Y.NRows; i++)
            //{
            //    RM += Math.Pow(Y[i, 0] - Summary.Mean(Y.GetColVector(0)), 2);
            //}
            //var F = (RM - rssUnResModel) / (bAux.NRows - 1) / (rssUnResModel / (Y.NRows - bAux.NRows));

            //return rssUnResModel;
            return temporal;
        }

        public static VarParameters VarTEURM(RMatrix X, int modelOrder)
        {
            var rss = 0.0;

            var xAux = new RMatrix(X);
            //var yAux = X.GetColVector(1);
            //for (int i = 0, j = 0; i < X.NCols; i++)
            //{
            //    if (i != 1)
            //    {
            //        //if (LeastSquaresTest(yAux, X.GetColVector(i)))
            //        //{

            //        //}
            //        xAux.InsertColVector(X.GetColVector(i), j);
            //        j++;
            //    }

            //}
            ////var modelOrder = ModelOrderGC(yIndex, uIndex, xAux)[1];

            ////xAux = RemoveColumn(yIndex, xAux);
            //xAux.InsertColVector(yAux, 0);
            var xActual = new List<double>(xAux.GetColVector(0).ToList());
            //var restrictedModel = new RMatrix();
            var unrestrictedModel = new RMatrix();

            for (var i = 0; i < xAux.NCols; i++)
            {
                //var xActual = new List<double>(xAux.GetColVector(i).ToList());
                var xPasada = new RMatrix(xAux.GetColVector(i).Length - modelOrder, modelOrder);
                var k = 0;

                for (var m = modelOrder - 1; m >= 0; m--)
                {
                    for (var j = 0; j < xAux.GetColVector(i).Length - modelOrder; j++)
                        xPasada[j, m] = xAux.GetColVector(i)[j + k];
                    k++;
                    //xActual.RemoveAt(0);
                    if (i == 0)
                    {
                        xActual.RemoveAt(0);
                    }
                }
                //Xt.InsertColVector(new RVector(xActual.ToArray()), 0);
                //xPasada.InsertColVector(new RVector(xActual.ToArray()), xPasada.NCols);
                for (var j = 0; j < xPasada.NCols; j++)
                {
                    unrestrictedModel.InsertColVector(xPasada.GetColVector(j), unrestrictedModel.NCols);
                    //if (i != uIndex)
                    //{
                    //    restrictedModel.InsertColVector(xPasada.GetColVector(j), restrictedModel.NCols);
                    //}
                }
            }
            var I = new RVector(unrestrictedModel.NRows);

            for (var i = 0; i < I.Length; i++)
            {
                I[i] = 1;
            }
            var Y = new RMatrix();
            Y.InsertColVector(new RVector(xActual.ToArray()), 0);
            unrestrictedModel.InsertColVector(I, 0);

            var bAux = RMatrix.Inverse(unrestrictedModel.GetTranspose() * unrestrictedModel) * unrestrictedModel.GetTranspose() * Y;

            var temporal = new RVector(unrestrictedModel.NRows);
            temporal = Transform(temporal, Y - unrestrictedModel * bAux);

            //var rss = 0.0;
            for (var i = 0; i < temporal.Length; i++)
            {
                rss += Math.Pow(temporal[i], 2);
            }
            //var cov = Summary.Covariance(Y - unrestrictedModel * bAux);
            var b = bAux.GetColVector(0).ToList();
            b.RemoveAt(0);

            var media = bAux[0, 0] / (1 - b.Sum());
            var sigma = 0.0;
            for (int i = modelOrder + 1; i < temporal.Length; i++)
            {
                sigma += Math.Pow(temporal[i], 2);
            }
            sigma = sigma / (temporal.Length - modelOrder);
            
            var varParam = new VarParameters(bAux[0,0], new RVector(b.ToArray()), temporal, sigma);

            unrestrictedModel = RemoveColumn(0, unrestrictedModel);
            for (int i = 0; i < unrestrictedModel.NCols; i++)
            {
                for (int j = 0; j < unrestrictedModel.NRows; j++)
                {
                    unrestrictedModel[j, i] = unrestrictedModel[j, i] - media;
                }
            }

            var covY = new RMatrix(unrestrictedModel.NCols, unrestrictedModel.NCols);
            for (int j = 0; j < unrestrictedModel.NCols; j++)
            {
                for (int l = 0; l < unrestrictedModel.NCols; l++)
                {
                    var vec = new RMatrix();
                    vec.InsertColVector(unrestrictedModel.GetColVector(j), 0);
                    vec.InsertColVector(unrestrictedModel.GetColVector(l), 1);
                    var temp = Summary.Covariance(vec);
                    if (j < l)
                        covY[j, l] = temp[0, 1];
                    else
                        covY[j, l] = temp[1, 0];
                }
            }
            var probY = 0.0;
            //for (int i = 0; i < unrestrictedModel.NRows; i++)
            //{
            //    var v = new RMatrix();
            //    v.InsertColVector(unrestrictedModel.GetRowVector(i), 0);

            //    var temp = v.GetTranspose() * RMatrix.Inverse(covY) * v;
            //    probY += Math.Pow(2 * Math.PI, -1.0 * modelOrder / 2)
            //        * Math.Pow(Math.Pow(sigma, -2), -1.0 * modelOrder / 2)
            //        * Math.Sqrt(RMatrix.QRDeterminant(RMatrix.Inverse(covY)))
            //        * Math.Exp(-0.5 * temp[0, 0] / Math.Pow(sigma, 2))
            //        ;

            //}
            for (int i = 0; i < temporal.Length; i++)
            {
                probY += Math.Exp(-0.5*Math.Pow(temporal[i]/sigma, 2))
                    /Math.Sqrt(2*Math.PI*Math.Pow(sigma, 2));
            }
            probY = probY / (Y.NRows);

            //var v = new RMatrix();
            //v.InsertColVector(unrestrictedModel.GetRowVector(0), 0);

            //var t = v.GetTranspose() * RMatrix.Inverse(covY) * v;
            //probY += Math.Pow(2 * Math.PI, -1.0 * modelOrder / 2)
            //    * Math.Pow(Math.Pow(sigma, -2), -1.0 * modelOrder / 2)
            //    * Math.Sqrt(RMatrix.QRDeterminant(RMatrix.Inverse(covY)))
            //    * Math.Exp(-0.5 * t[0, 0] / Math.Pow(sigma, 2))
            //    ;
            //for (var i = modelOrder; i < temporal.Length; i++)
            //{
            //    probY *= Math.Exp(-0.5 * Math.Pow(temporal[i] / sigma, 2))
            //        / Math.Sqrt(2 * Math.PI * Math.Pow(sigma, 2));
            //}

            //return rss;
            //return RMatrix.QRDeterminant(cov);
            return new VarParameters(bAux[0, 0], new RVector(b.ToArray()), temporal, sigma, unrestrictedModel);
        }

        public static VarParameters JointPDF(RMatrix X, int modelOrder)
        {
            var rss = 0.0;
            var xAux = new RMatrix(X);
            var xActual = new List<double>(xAux.GetColVector(0).ToList());
            var unrestrictedModel = new RMatrix();

            for (var i = 0; i < xAux.NCols; i++)
            {
                var xPasada = new RMatrix(xAux.GetColVector(i).Length - modelOrder, modelOrder);
                var k = 0;

                for (var m = modelOrder - 1; m >= 0; m--)
                {
                    for (var j = 0; j < xAux.GetColVector(i).Length - modelOrder; j++)
                        xPasada[j, m] = xAux.GetColVector(i)[j + k];
                    k++;
                    if (i == 0)
                    {
                        xActual.RemoveAt(0);
                    }
                }
                for (var j = 0; j < xPasada.NCols; j++)
                {
                    unrestrictedModel.InsertColVector(xPasada.GetColVector(j), unrestrictedModel.NCols);

                }
            }
            var I = new RVector(unrestrictedModel.NRows);

            for (var i = 0; i < I.Length; i++)
            {
                I[i] = 1;
            }
            var Y = new RMatrix();
            Y.InsertColVector(new RVector(xActual.ToArray()), 0);
            unrestrictedModel.InsertColVector(I, 0);

            var bAux = RMatrix.Inverse(unrestrictedModel.GetTranspose() * unrestrictedModel) * unrestrictedModel.GetTranspose() * Y;

            var temporal = new RVector(unrestrictedModel.NRows);
            temporal = Transform(temporal, Y - unrestrictedModel * bAux);

            for (var i = 0; i < temporal.Length; i++)
            {
                rss += Math.Pow(temporal[i], 2);
            }
            var b = bAux.GetColVector(0).ToList();
            b.RemoveAt(0);

            var media = bAux[0, 0] / (1 - b.Sum());
            var sigma = 0.0;
            for (int i = modelOrder + 1; i < temporal.Length; i++)
            {
                sigma += Math.Pow(temporal[i], 2);
            }
            sigma = sigma / (temporal.Length - modelOrder);

            unrestrictedModel = RemoveColumn(0, unrestrictedModel);
            for (int i = 0; i < unrestrictedModel.NCols; i++)
            {
                for (int j = 0; j < unrestrictedModel.NRows; j++)
                {
                    unrestrictedModel[j, i] = unrestrictedModel[j, i] - media;
                }
            }

            var covY = new RMatrix(unrestrictedModel.NCols, unrestrictedModel.NCols);
            for (int j = 0; j < unrestrictedModel.NCols; j++)
            {
                for (int l = 0; l < unrestrictedModel.NCols; l++)
                {
                    var vec = new RMatrix();
                    vec.InsertColVector(unrestrictedModel.GetColVector(j), 0);
                    vec.InsertColVector(unrestrictedModel.GetColVector(l), 1);
                    var temp = Summary.Covariance(vec);
                    if (j < l)
                        covY[j, l] = temp[0, 1];
                    else
                        covY[j, l] = temp[1, 0];
                }
            }
            var probY = 0.0;
            for (int i = 0; i < unrestrictedModel.NRows; i++)
            {
                var v = new RMatrix();
                v.InsertColVector(unrestrictedModel.GetRowVector(i), 0);

                var temp = v.GetTranspose() * RMatrix.Inverse(covY) * v;
                probY += Math.Pow(2 * Math.PI, -1.0 * modelOrder / 2)
                    * Math.Pow(Math.Pow(sigma, -2), -1.0 * modelOrder / 2)
                    * Math.Sqrt(RMatrix.QRDeterminant(RMatrix.Inverse(covY)))
                    * Math.Exp(-0.5 * temp[0, 0] / Math.Pow(sigma, 2))
                    ;

            }

            probY = probY / (Y.NRows * sigma);

            return new VarParameters(bAux[0, 0], new RVector(b.ToArray()), temporal, sigma, unrestrictedModel);
        }

        public static VarParameters VarTERM(IReadOnlyList<double> y, int modelOrder)
        {
            //if (u == null) throw new ArgumentNullException(nameof(u));
            if (y == null) throw new ArgumentNullException(nameof(y));
            
            var rss = 0.0;
            //var aic = new double[20];
            var N = 0;

            //var histogram = new RMatrix(Summary.Histogram(new RVector(y.ToArray()), (int)Math.Sqrt(y.Count)));
            //var Y = new RVector(histogram.NCols);
            //for (var i = 0; i < histogram.NCols; i++)
            //{
            //    Y[i] = histogram[0, i];
            //}
            var xActual = new List<double>(y);
            var xPasada = new RMatrix(y.Count - modelOrder, modelOrder);
            var Y = new RMatrix();

            //var temporal = new RVector();
            var k = 0;

            for (var i = modelOrder - 1; i >= 0; i--)
            {
                for (var j = 0; j < y.Count - modelOrder; j++)
                    xPasada[j, i] = y[j + k];
                k++;
                xActual.RemoveAt(0);
            }
            Y.InsertColVector(new RVector(xActual.ToArray()), 0);

            //var b = OrdinalLeastSquares(Y, xPasada);
            //var bAux = new RMatrix(xPasada.NCols + 1, 1);
            var I = new RVector(xPasada.NRows);
            var temporal = new RVector(xPasada.NRows);
            for (int i = 0; i < I.Length; i++)
            {
                I[i] = 1;
            }
            xPasada.InsertColVector(I, 0);
            var bAux = RMatrix.Inverse(xPasada.GetTranspose() * xPasada) * xPasada.GetTranspose() * Y;

            //RSS = Transform(RSS, Y - xPasada * bAux);

            temporal = Transform(temporal, Y - xPasada * bAux);
            for (var m = 0; m < temporal.Length; m++)
            {
                rss += Math.Pow(temporal[m], 2);
            }

            var b = bAux.GetColVector(0).ToList();
            b.RemoveAt(0);

            var media = bAux[0, 0]/(1 - b.Sum());
            var sigma = 0.0;
            for (int i = modelOrder + 1; i < temporal.Length; i++)
            {
                sigma += Math.Pow(temporal[i], 2);
            }
            sigma = sigma/(temporal.Length - modelOrder);
            

            //var prob = Gauss(sigma/(1 - Math.Pow(b[0], 2)), temporal[0]);
            //var Likelihood = Math.Log(prob, 2);
            //for (var i = 1; i < temporal.Length; i++)
            //{
            //    Likelihood += Math.Log(Gauss(sigma, temporal[i]), 2);
            //}
            xPasada = RemoveColumn(0, xPasada);
            for (int i = 0; i < xPasada.NCols; i++)
            {
                for (int j = 0; j < xPasada.NRows; j++)
                {
                    xPasada[j, i] = xPasada[j, i] - media;
                }
            }

            var covY = new RMatrix(xPasada.NCols, xPasada.NCols);
            for (int j = 0; j < xPasada.NCols; j++)
            {
                for (int l = 0; l < xPasada.NCols; l++)
                {
                    var vec = new RMatrix();
                    vec.InsertColVector(xPasada.GetColVector(j), 0);
                    vec.InsertColVector(xPasada.GetColVector(l), 1);
                    var temp = Summary.Covariance(vec);
                    if (j < l)
                        covY[j, l] = temp[0, 1];
                    else
                        covY[j, l] = temp[1, 0];
                }
            }
            var probY = 0.0;
            //for (var i = 0; i < xPasada.NRows; i++)
            //{
            //    var v = new RMatrix();
            //    v.InsertColVector(xPasada.GetRowVector(i), 0);

            //    var temp = v.GetTranspose() * RMatrix.Inverse(covY) * v;
            //    probY += Math.Pow(2 * Math.PI, -1.0 * modelOrder / 2)
            //        * Math.Pow(Math.Pow(sigma, -2), -1.0 * modelOrder / 2)
            //        * Math.Sqrt(RMatrix.QRDeterminant(RMatrix.Inverse(covY)))
            //        * Math.Exp(-0.5 * temp[0, 0] / Math.Pow(sigma, 2))
            //        ;

            //}

            for (int i = 0; i < temporal.Length; i++)
            {
                probY += Math.Exp(-0.5 * Math.Pow(temporal[i] / sigma, 2))
                    / Math.Sqrt(2 * Math.PI * Math.Pow(sigma, 2));
            }

            probY = probY / (Y.NRows);

            //var v = new RMatrix();
            //v.InsertColVector(xPasada.GetRowVector(0), 0);

            //var t = v.GetTranspose() * RMatrix.Inverse(covY) * v;
            //probY += Math.Pow(2 * Math.PI, -1.0 * modelOrder / 2)
            //    * Math.Pow(Math.Pow(sigma, -2), -1.0 * modelOrder / 2)
            //    * Math.Sqrt(RMatrix.QRDeterminant(RMatrix.Inverse(covY)))
            //    * Math.Exp(-0.5 * t[0, 0] / Math.Pow(sigma, 2))
            //    ;
            //for (var i = modelOrder; i < temporal.Length; i++)
            //{
            //    probY *= Math.Exp(-0.5 * Math.Pow(temporal[i] / sigma, 2))
            //        / Math.Sqrt(2 * Math.PI * Math.Pow(sigma, 2));
            //}

            //for (var i = modelOrder; i < temporal.Length; i++)
            //{
            //    probY -= 0.5 * Math.Pow(temporal[i] / sigma, 2);
            //}
            //var Likelihood = -0.5 * (temporal.Length - modelOrder) * Math.Log(2 * Math.PI)
            //        - 0.5 * (temporal.Length - modelOrder) * Math.Log(Math.Pow(sigma, 2))
            //        - probY;


            //var Likelihood = Math.Log(probY);
            //var cov = Summary.Covariance(Y - xPasada*bAux);
            //aic[n - 1] = Math.Log(rss[n - 1]) + (double)2 * n / Y.NRows;
            //for (var i = 1; i < aic.Length; i++)
            //{
            //    if (aic[i] < aic[N])
            //    {
            //        N = i;
            //    }
            //}
            //N = N + 1;

            //return N;
            //return rss;
            //return RMatrix.QRDeterminant(cov);


            return new VarParameters(bAux[0, 0], new RVector(b.ToArray()), temporal, sigma, xPasada);
        }

        public static RMatrix LinearRegression(RMatrix X, int[,] retardo)
        {
            var B = new RMatrix();
            var Bi = new RVector(X.NCols);
            var yDis = new List<List<double>>();
            for (var i = 0; i < X.NCols; i++)
            {
                var Y = new RMatrix();
                var xAux = new RMatrix(X);  
                Y.InsertColVector(X.GetColVector(i), 0);
                xAux = RemoveColumn(i, xAux);
                Bi = Transform(Bi, OrdinalLeastSquares(Y, xAux));
                
                B.InsertRowVector(Bi, i);
                //var h = Discretization(Y, xAux, Bi);
                //yDis.Add(h.ToList());
                
            }
            //var c = TransferEntropy(yDis, retardo);
            return B;
        }

        public static RMatrix VectorAutoRegession(RMatrix X)
        {
            var B = new RMatrix();
            var Bi = new RVector(X.NCols + 1);
            var yDis = new List<List<double>>();
            
            for (var i = 0; i < X.NCols; i++)
            {
                var Y = new RMatrix();
                var xAux = new RMatrix(X);
                Y.InsertColVector(X.GetColVector(i), 0);
                xAux = RemoveColumn(i, xAux);
                var Yt = new RMatrix(Y);
                Y = RemoveRow(0, Y);
                xAux.InsertColVector(Yt.GetColVector(0), 0);
                xAux = RemoveRow(xAux.NRows - 1, xAux);
                Bi = Transform(Bi, OrdinalLeastSquares(Y, xAux));

                B.InsertRowVector(Bi, i);
                //var h = Discretization(Y, xAux, Bi);
                //yDis.Add(h.ToList());

            }
            //var c = TransferEntropy(yDis, retardo);
            return B;
        }

        public static RMatrix VAR(RMatrix X, int modelOrder)
        {
            var B = new RMatrix();
            for (var m = 0; m < X.NCols; m++)
            {
                var y = X.GetColVector(m).ToList();

                var xActual = new List<double>(y);
                var xPasada = new RMatrix(y.Count - modelOrder, modelOrder);
                var Xt = new RMatrix(X);
                Xt = RemoveColumn(m, Xt);
                var Y = new RMatrix();
                var k = 0;

                for (var i = modelOrder - 1; i >= 0; i--)
                {
                    for (var j = 0; j < y.Count - modelOrder; j++)
                        xPasada[j, i] = y[j + k];
                    k++;
                    xActual.RemoveAt(0);
                    Xt = RemoveRow(Xt.NRows - 1, Xt);
                }
                Y.InsertColVector(new RVector(xActual.ToArray()), 0);
                
                for (var i = 0; i < Xt.NCols; i++)
                {
                    xPasada.InsertColVector(Xt.GetColVector(i), xPasada.NCols);

                }

                
                var Bi = OrdinalLeastSquares(Y, xPasada);

                B.InsertRowVector(Bi.GetColVector(0), index: m);
            }
            return B;
        }

        public static RMatrix DirectTransferEntropy(RMatrix Discret, RMatrix TE)
        {
            var CausalMatrixDirect = new double[Discret.NCols, Discret.NCols];
            var CausalMatrixIndirect = new double[Discret.NCols, Discret.NCols];
            var FinalCausalMatrix = new double[Discret.NCols, Discret.NCols];
            var ETE = new double[Discret.NCols, Discret.NCols];
            for (int i = 0; i < Discret.NCols; i++)
            {
                for (int j = 0; j < Discret.NCols; j++)
                {
                    for (int k = 0; k < Discret.NCols; k++)
                    {
                        if (i != j && i != k && j!=k)
                        {
                            if (!Equals(TE[i, k], 0.0) && !Equals(TE[k, j], 0.0))
                            {
                                var tempZY =
                                    DTExyz(u: Discret.GetColVector(k).ToList(), y: Discret.GetColVector(j).ToList(),
                                        z: Discret.GetColVector(i).ToList());

                                if (!Equals(tempZY.Coeficiente, 0.0))
                                {
                                    var tempXY = DTExyz(Discret.GetColVector(i).ToList(), Discret.GetColVector(j).ToList(),
                                        Discret.GetColVector(k).ToList());

                                    CausalMatrixDirect[k, j] = tempZY.Coeficiente;
                                    CausalMatrixIndirect[i, j] = tempXY.Coeficiente;

                                }
                            }
                        }
                    }
                }
            }

            //for (int i = 0; i < Discret.NCols; i++)
            //{
            //    for (int j = 0; j < Discret.NCols; j++)
            //    {
            //        FinalCausalMatrix[i, j] = CausalMatrixDirect[i, j] > CausalMatrixIndirect[i,j]
            //            ? CausalMatrixDirect[i, j]
            //            : CausalMatrixIndirect[i, j];
            //    }
            //}
            for (int i = 0; i < Discret.NCols; i++)
            {
                for (int j = 0; j < Discret.NCols; j++)
                {
                    ETE[i, j] = CausalMatrixIndirect[i, j] - CausalMatrixIndirect[j, i];
                    if (ETE[i, j] < 0.001)
                    {
                        ETE[i, j] = 0;
                    }
                }
            }

            return new RMatrix(ETE);
        }

        public static RMatrix TransferEntropy(RMatrix Discret)
        {
            var CoefTE = new double[Discret.NCols, Discret.NCols];
            var ETE = new double[Discret.NCols, Discret.NCols];

            //Discret.ZScoreNormalize();
            //var TE = new double[10];
            //var n1 = new double[10];
            //var n2 = new double[10];
            for (var i = 0; i < Discret.NCols; i++)
            {
                for (int j = 0; j < Discret.NCols; j++)
                {
                    if (i != j)
                    {
                        //var X = new RMatrix();
                        //X.InsertColVector(Discret.GetColVector(j), 0);
                        //X.InsertColVector(Discret.GetColVector(i), 1);
                        //var ERM = VarTERM(X.GetColVector(0).ToList(), 4);
                        //var EURM = VarTEURM(X, 4);
                        //if (F_Test(ERM, EURM, 4, X.GetColVector(0).Length, X.NCols))
                        //    CoefTE[i, j].Coeficiente = 0.5*Math.Log(ERM/EURM);
                        //else
                        //    CoefTE[i, j].Coeficiente = 0;

                        //for (var k = 1; k < 11; k++)
                        //{

                        //    var ERM = VarTERM(X.GetColVector(0).ToList(), k);
                        //    var EURM = VarTEURM(X, k);

                        //    n1[k - 1] = Math.Log(ERM) + 2 * k * 1 / X.NRows;
                        //    n2[k - 1] = Math.Log(EURM) + 2 * k * 2 / X.NRows;

                        //    if (F_Test(ERM, EURM, k, X.GetColVector(0).Length, X.NCols))
                        //        TE[k - 1] = (1 - EURM / ERM);
                        //    else
                        //        TE[k - 1] = 0;
                        //}

                        //var n = ModelOrder(Discret.GetColVector(i).ToList(), Discret.GetColVector(j).ToList());

                        //var X = new RMatrix();
                        //X.InsertColVector(Discret.GetColVector(j), 0);
                        //X.InsertColVector(Discret.GetColVector(i), 1);
                        //var ERM = new double[10];
                        //var EURM = new double[10];
                        //for (var k = 1; k < 11; k++)
                        //{
                        //    ERM[k - 1] = VarTERM(X.GetColVector(0).ToList(), k);
                        //    EURM[k - 1] = VarTEURM(X, k);
                        //}
                        //var n = Aikake(ERM, EURM, X.GetColVector(0).Length, X.NCols);
                        //if (F_Test(ERM[n[1] - 1], EURM[n[1] - 1], n[1], X.GetColVector(0).Length, X.NCols))
                        //{
                        //    CoefTE[i, j].Coeficiente = (1 - EURM[n[1] - 1] / ERM[n[1] - 1]);
                        //    //CoefTE[i, j].Coeficiente = CoefTE[i, j].Coeficiente / Hy(X.GetColVector(0));
                        //    CoefTE[i, j].Lambda = n[0];
                        //}
                        //else
                        //{
                        //    CoefTE[i, j].Coeficiente = 0;
                        //    CoefTE[i, j].Lambda = n[1];
                        //}

                        var temp = NTE(Discret.GetColVector(i).ToList(), Discret.GetColVector(j).ToList());
                        //if (temp.Coeficiente < 0.01)
                        //    continue;

                        CoefTE[i, j] = temp.Coeficiente;
                        //if (CoefTE[i, j].coeficiente < (1.85 * Math.Pow(Discret.NRows, -0.41) + (2.37 * Math.Pow(Discret.NRows, -0.53))))
                        //{
                        //    CoefTE[i, j].coeficiente = 0;
                        //}

                    }

                }
            }
            for (int i = 0; i < Discret.NCols; i++)
            {
                for (int j = 0; j < Discret.NCols; j++)
                {
                    ETE[i, j] = CoefTE[i, j] - CoefTE[j, i];
                    if (ETE[i, j] > 0.001) continue;
                    ETE[i, j] = 0;
                }
            }
            return new RMatrix(CoefTE);
        }

        public static int Signo(double a)
        {
            if (a < 0)
                return -1;
            return 1;
        }

        public static RVector Cuntify(IReadOnlyList<double> y, List<double> codeBook)
        {
            var Y = new RVector(y.Count);
            for (var k = 0; k < y.Count; k++)
            {
                for (var i = 0; i < codeBook.Count - 1; i++)
                {
                    if (y[k] >= codeBook[i] && y[k] <= codeBook[i + 1])
                    {
                        Y[k] = codeBook[i];
                    }
                }
            }

            return Y;
        }

        public static RVector Cuntify(IReadOnlyList<double> y, int cantidad)
        {
            //var count = (int)Math.Sqrt(y.Count)/3;
            var Y = new RVector(y.Count);
            var min = y.Min();
            var max = y.Max();
            var len = Math.Abs(max) + Math.Abs(min);
            var codeBook = new RVector(cantidad);
            var L = new RVector(cantidad + 1);
            var h = len/(L.Length - 1);
            L[0] = min;
            for (var i = 1; i < L.Length; i++)
            {
                L[i] = L[i - 1] + h;
                codeBook[i - 1] = L[i - 1] + h/2;
            }
            for (var k = 0; k < Y.Length; k++)
            {
                for (var i = 1; i < L.Length; i++)
                {
                    if (y[k] >= L[i - 1] && y[k] <= L[i])
                    {
                        Y[k] = codeBook[i - 1];
                    }
                }
            }


            //for (int i = 0, u = 0; i < count; i++)
            //{
            //    Y[i] = Math.Round((count - 1) * (y[u] - y.Min()) / (y.Max() - y.Min())) + 1;
            //    u += y.Count / count;
            //}

            return Y;

            //var Mu = 255;
            //var Y = new List<double>();
            //foreach (var t in y)
            //{
            //    Y.Add(y.Max()
            //          * (Math.Log(1 + Mu * Math.Abs(t) / y.Max()) / Math.Log(1 + Mu))
            //          * Signo(t))
            //        ;
            //}


            //return new RVector(Y.ToArray());
        }

        public static List<double> LloydAlgorithm(IReadOnlyList<double> u, IReadOnlyList<double> y)
        {
            var finalcodeBook = new List<double>();
            var epsilon = 0.000001;
            var Y = Cuntify(y, 4);
            //var Y = new RVector(y.ToArray());
            var tempCodeBook = new List<double> { Y[0] };
            
            foreach (var t in Y.ToList().Where(t => !tempCodeBook.Contains(t)))
                tempCodeBook.Add(t);
            tempCodeBook.Sort();

            var lista = new List<List<double>>();
            for (var i = 0; i < tempCodeBook.Count - 1; i++)
            {
                var temp = new List<double>();
                foreach (var t in y)
                    if (Math.Sqrt(Math.Pow(t - tempCodeBook[i], 2))
                        < Math.Sqrt(Math.Pow(tempCodeBook[i] - tempCodeBook[i + 1], 2)))
                        temp.Add(t);
                lista.Add(temp);
            }

            var freqY = new List<double>(tempCodeBook.Count);
            freqY.AddRange(tempCodeBook.Select(t1 => Y.ToList().Count(t => Equals(t1, t))).Select(count => (double)count));

            //var mediaY = y.Average();
            //var sigmaY = Summary.StandardDeviation(new RVector(y.ToArray()));
            var D = 0.0;
            var D1 = 0.0;
            foreach (var t in y)
            {
                var index = 0;
                for (var i = 0; i < lista.Count; i++)
                    if (lista[i].Contains(t))
                        index = i;
                var c = Math.Pow(t - tempCodeBook[index], 2);
                //var probY = Math.Exp(-0.5 * Math.Pow(t - mediaY, 2) / Math.Pow(sigmaY, 2))
                //            / (Math.Sqrt(2 * Math.PI) * sigmaY);
                var probY = freqY[index]/y.Count;

                D += probY * c;
            }

            do
            {
                var newCodeBook = new List<double>();
                for (var i = 0; i < tempCodeBook.Count - 1; i++)
                {
                    var L = tempCodeBook[i];
                    var R = tempCodeBook[i + 1];
                    var index = 0;
                    var py = freqY[i] / y.Count;
                    double codeWord = 0;
                    for (int j = 0; j < y.Count; j++)
                    {
                        if (Y[j] >= L && Y[j] <= R) codeWord += y[j] * py;
                    }
                    newCodeBook.Add(Math.Round(codeWord));
                }
                newCodeBook.Sort();

                var lista1 = new List<List<double>>();
                for (var i = 0; i < newCodeBook.Count - 1; i++)
                {
                    var temp = new List<double>();
                    foreach (var t in y)
                        if (Math.Sqrt(Math.Pow(t - newCodeBook[i], 2))
                            < Math.Sqrt(Math.Pow(newCodeBook[i] - newCodeBook[i + 1], 2)))
                            temp.Add(t);
                    lista1.Add(temp);
                }

                foreach (var t in y)
                {
                    var index = 0;
                    for (var i = 0; i < lista1.Count; i++)
                        if (lista1[i].Contains(t))
                            index = i;
                    var c = Math.Pow(t - newCodeBook[index], 2);
                    var probY = freqY[index] / y.Count;

                    D1 += probY * c;
                }
            } while (Math.Abs((D - D1)/D) < epsilon);

            return finalcodeBook;
        }

        public static RVector RemoveAt(int pos, RVector vector)
        {
            var temp = vector.ToList();
            temp.RemoveAt(pos);
            return new RVector(temp.ToArray());
        }

        public static Correlacion DTExyz(IReadOnlyList<double> u, IReadOnlyList<double> y, IReadOnlyList<double> z)
        {
            if (u == null) throw new ArgumentNullException(nameof(u));
            if (y == null) throw new ArgumentNullException(nameof(y));
            if (z == null) throw new ArgumentNullException(nameof(z));

            var coef = new Correlacion { Coeficiente = 0 };
            var len = 4;
            var Y = Cuntify(y, len).ToList();
            var U = Cuntify(u, len).ToList();
            var Z = Cuntify(z, len).ToList();

            //var codeBook = new List<double>();
            //var min = y.Min();
            //var max = y.Max();
            //var h = (max - min) / 5;
            //codeBook.Add(min);
            //for (int i = 0; i < 5; i++)
            //{
            //    codeBook.Add(Math.Round(codeBook[i] + h, 2));
            //}

            //var Y = Cuntify(y, codeBook).ToList();
            //var U = Cuntify(u, codeBook).ToList();
            //var Z = Cuntify(z, codeBook).ToList();

            //var n = UnivarModelOrder(y);
            //var n1 = n + 1;

            //var X = new RMatrix();
            //X.InsertColVector(Y, 0);
            //X.InsertColVector(U, 1);
            //var l = MultiModelOrder(0, X);

            //X = RemoveColumn(1, X);
            //X.InsertColVector(Z, 1);
            //var k = MultiModelOrder(0, X);

            //var l = UnivarModelOrder(u);
            //var k = UnivarModelOrder(z);
            var n = 1;
            var n1 = n + 1;
            var l = n;
            var k = n;

            var DTE = new double[30];

            var codeBook = new List<double> { Y[0] };

            foreach (var t in Y.ToList().Where(t => !codeBook.Contains(t)))
                codeBook.Add(t);

            foreach (var t in U.ToList().Where(t => !codeBook.Contains(t)))
                codeBook.Add(t);

            foreach (var t in Z.ToList().Where(t => !codeBook.Contains(t)))
                codeBook.Add(t);

            codeBook.Sort();

            var freqY = new List<int>(codeBook.Count);
            var freqU = new List<int>(codeBook.Count);
            var freqZ = new List<int>(codeBook.Count);

            freqY.AddRange(codeBook.Select(t1 => Y.Count(t => Equals(t1, t))).Select(count => count));

            freqU.AddRange(codeBook.Select(t1 => U.Count(t => Equals(t1, t))).Select(count => count));

            freqZ.AddRange(codeBook.Select(t1 => Z.Count(t => Equals(t1, t))).Select(count => count));

            var Hy = freqY.Where(t => !Equals(t, 0)).Aggregate(0.0, (current, t) => current - (double)t / Y.Count * Math.Log((double)t / Y.Count, 2));
            var Hu = freqU.Where(t => !Equals(t, 0)).Aggregate(0.0, (current, t) => current - (double)t / U.Count * Math.Log((double)t / U.Count, 2));
            var Hz = freqU.Where(t => !Equals(t, 0)).Aggregate(0.0, (current, t) => current - (double)t / Z.Count * Math.Log((double)t / U.Count, 2));

            for (var lamda = 1; lamda < 31; lamda++)
            {
                
                var freqYn1 = new List<int>();
                var freqYn1UlZk = new List<int>();
                var freqYnUlZk = new List<int>();
                var freqYnZk = new List<int>();
                var freqYn1Zk = new List<int>();
                var freqYn1Ul = new List<int>();
                var freqYnUl = new List<int>();
                var freqYn = new List<int>();
                var freqUn = new List<int>();
                /*********************************************/
                foreach (var ty in codeBook)
                {
                    foreach (var tz in codeBook)
                    {
                        foreach (var tu in codeBook)
                        {
                            var cdY = new double[n];
                            var cdU = new double[l];
                            var cdZ = new double[k];
                            for (var i = 0; i < cdY.Length; i++)
                                cdY[i] = ty;
                            for (var i = 0; i < cdU.Length; i++)
                                cdU[i] = tu;
                            for (var i = 0; i < cdZ.Length; i++)
                                cdZ[i] = tz;

                            var count = 0;

                            for (int i = 0, m = 0, d = 0; (i < Y.Count - n && m < U.Count - l && d < Z.Count - k);)
                            {
                                var tempY = new double[n];
                                for (var j = 0; j < n; j++)
                                    tempY[j] = Y[i + j];

                                var tempU = new double[l];
                                for (var j = 0; j < l; j++)
                                    tempU[j] = U[i + j];

                                var tempZ = new double[k];
                                for (var j = 0; j < k; j++)
                                    tempZ[j] = Z[i + j];

                                if (Equals(new RVector(tempY), new RVector(cdY))
                                    && Equals(new RVector(tempU), new RVector(cdU))
                                    && Equals(new RVector(tempZ), new RVector(cdZ)))
                                {
                                    count++;
                                }
                                i++;
                                m++;
                                d++;
                            }

                            freqYnUlZk.Add(count);
                        }
                    }
                }
                /*********************************************/
                foreach (var ty in codeBook)
                {
                    foreach (var tz in codeBook)
                    {
                        foreach (var tu in codeBook)
                        {
                            var cdY = new double[n1];
                            var cdU = new double[l];
                            var cdZ = new double[k];
                            for (var i = 0; i < cdY.Length; i++)
                                cdY[i] = ty;
                            for (var i = 0; i < cdU.Length; i++)
                                cdU[i] = tu;
                            for (var i = 0; i < cdZ.Length; i++)
                                cdZ[i] = tz;

                            var count = 0;

                            for (int i = 0, m = 0, d = 0; (i < Y.Count - n1 && m < U.Count - l && d < Z.Count - k);)
                            {
                                var tempY = new double[n1];
                                for (var j = 0; j < n1; j++)
                                    tempY[j] = Y[i + j];

                                var tempU = new double[l];
                                for (var j = 0; j < l; j++)
                                    tempU[j] = U[i + j];

                                var tempZ = new double[k];
                                for (var j = 0; j < k; j++)
                                    tempZ[j] = Z[i + j];

                                if (Equals(new RVector(tempY), new RVector(cdY))
                                    && Equals(new RVector(tempU), new RVector(cdU))
                                    && Equals(new RVector(tempZ), new RVector(cdZ)))
                                {
                                    count++;
                                }
                                i++;
                                m++;
                                d++;
                            }

                            freqYn1UlZk.Add(count);
                        }
                    }
                }
                /*********************************************/
                foreach (var ty in codeBook)
                {
                    foreach (var tz in codeBook)
                    {
                        foreach (var tu in codeBook)
                        {
                            var cdY = new double[n];
                            var cdZ = new double[k];
                            for (var i = 0; i < cdY.Length; i++)
                                cdY[i] = ty;
                            for (var i = 0; i < cdZ.Length; i++)
                                cdZ[i] = tz;

                            var count = 0;

                            for (int i = 0, m = 0; (i < Y.Count - n && m < U.Count - k);)
                            {
                                var tempY = new double[n];
                                for (var j = 0; j < n; j++)
                                    tempY[j] = Y[i + j];

                                var tempU = new double[k];
                                for (var j = 0; j < k; j++)
                                    tempU[j] = U[i + j];

                                if (Equals(new RVector(tempY), new RVector(cdY))
                                    && Equals(new RVector(tempU), new RVector(cdZ)))
                                {
                                    count++;
                                }
                                i++;
                                m++;
                            }

                            freqYnZk.Add(count);
                        }
                    }
                }
                /*********************************************/
                foreach (var ty in codeBook)
                {
                    foreach (var tz in codeBook)
                    {
                        foreach (var tu in codeBook)
                        {
                            var cdY = new double[n1];
                            var cdZ = new double[k];
                            for (var i = 0; i < cdY.Length; i++)
                                cdY[i] = ty;
                            for (var i = 0; i < cdZ.Length; i++)
                                cdZ[i] = tz;

                            var count = 0;

                            for (int i = 0, m = 0; (i < Y.Count - n1 && m < U.Count - k);)
                            {
                                var tempY = new double[n1];
                                for (var j = 0; j < n1; j++)
                                    tempY[j] = Y[i + j];

                                var tempU = new double[k];
                                for (var j = 0; j < k; j++)
                                    tempU[j] = U[i + j];

                                if (Equals(new RVector(tempY), new RVector(cdY))
                                    && Equals(new RVector(tempU), new RVector(cdZ)))
                                {
                                    count++;
                                }
                                i++;
                                m++;
                            }

                            freqYn1Zk.Add(count);
                        }
                    }
                }
                /*********************************************/
                //foreach (var t in codeBook)
                //{
                //    var cdY = new double[n];
                //    var cdU = new double[l];
                //    for (var i = 0; i < cdY.Length; i++)
                //        cdY[i] = t;
                //    for (var i = 0; i < cdU.Length; i++)
                //        cdU[i] = t;

                //    var count = 0;

                //    for (int i = 0, m = 0; (i < Y.Length / n && m < U.Length / l);)
                //    {
                //        var tempY = new double[n];
                //        for (var j = 0; j < n; j++)
                //            tempY[j] = Y[i + j];

                //        var tempU = new double[l];
                //        for (var j = 0; j < l; j++)
                //            tempU[j] = U[i + j];

                //        if (Equals(new RVector(tempY), new RVector(cdY))
                //            && Equals(new RVector(tempU), new RVector(cdU)))
                //        {
                //            count++;
                //        }
                //        i += n;
                //        m += l;
                //    }

                //    freqYnUl.Add(count);
                //}
                /*********************************************/
                //foreach (var t in codeBook)
                //{
                //    var cdYn1 = new double[n1];
                //    var cdU = new double[l];
                //    for (var i = 0; i < cdYn1.Length; i++)
                //        cdYn1[i] = t;
                //    for (var i = 0; i < cdU.Length; i++)
                //        cdU[i] = t;

                //    var count = 0;

                //    for (int i = 0, m = 0; (i < Y.Length / n1 && m < U.Length / l);)
                //    {
                //        var tempYn1 = new double[n1];
                //        for (var j = 0; j < n1; j++)
                //            tempYn1[j] = Y[i + j];

                //        var tempU = new double[l];
                //        for (var j = 0; j < l; j++)
                //            tempU[j] = U[i + j];

                //        if (Equals(new RVector(tempYn1), new RVector(cdYn1))
                //            && Equals(new RVector(tempU), new RVector(cdU)))
                //        {
                //            count++;
                //        }
                //        i += n1;
                //        m += l;
                //    }

                //    freqYn1Ul.Add(count);
                //}
                /*********************************************/
                foreach (var ty in codeBook)
                {
                    foreach (var tz in codeBook)
                    {
                        foreach (var tu in codeBook)
                        {
                            var cdYn1 = new double[n1];
                            for (var i = 0; i < cdYn1.Length; i++)
                            {
                                cdYn1[i] = ty;
                            }
                            var count = 0;

                            for (int i = 0; i < Y.Count - n1;)
                            {
                                var temp = new double[n1];
                                for (var j = 0; j < n1; j++)
                                {
                                    temp[j] = Y[i + j];
                                }
                                if (Equals(new RVector(temp), new RVector(cdYn1)))
                                {
                                    count++;
                                }
                                i++;
                            }
                            freqYn1.Add(count);
                        }
                    }
                }
                /*********************************************/
                foreach (var ty in codeBook)
                {
                    foreach (var tz in codeBook)
                    {
                        foreach (var tu in codeBook)
                        {
                            var cdYn = new double[n];
                            for (var i = 0; i < cdYn.Length; i++)
                            {
                                cdYn[i] = ty;
                            }
                            var count = 0;

                            for (int i = 0; i < Y.Count - n;)
                            {
                                var temp = new double[n];
                                for (var j = 0; j < n; j++)
                                {
                                    temp[j] = Y[i + j];
                                }
                                if (Equals(new RVector(temp), new RVector(cdYn)))
                                {
                                    count++;
                                }
                                i++;
                            }
                            freqYn.Add(count);
                        }
                    }
                }
                /*********************************************/

                //var Hyn = freqYn.Where(t => !Equals(t, 0)).Aggregate(0.0, (current, t) => current - (double)t / Y.ToList().Count * Math.Log((double)t / Y.ToList().Count, 2));
                //var Hyn1 = freqYn1.Where(t => !Equals(t, 0)).Aggregate(0.0, (current, t) => current - (double)t / Y.ToList().Count * Math.Log((double)t / Y.ToList().Count, 2));
                var Hyn1 = 0.0;
                //var Hynulzk = freqYnUlZk.Where(t => !Equals(t, 0)).Aggregate(0.0, (current, t) => current - (double)t / Y.ToList().Count * Math.Log((double)t / Y.ToList().Count, 2));
                //var Hyn1ulzk = freqYn1UlZk.Where(t => !Equals(t, 0)).Aggregate(0.0, (current, t) => current - (double)t / Y.ToList().Count * Math.Log((double)t / Y.ToList().Count, 2));
                var Hyn1ulzk = 0.0;
                /*********************************************/
                //var probJoint = 0.0;
                var probYnZk = 0.0;
                var probYn1Zk = 0.0;
                //var probY = 0.0;
                //var probU = 0.0;
                var probYn1 = 0.0;
                var probYn = 0.0;
                var probYnUlZk = 0.0;
                var probYn1UlZk = 0.0;
                //var probYnUl = 0.0;
                //var probYn1Ul = 0.0;
                /*********************************************/
                for (int i = 0; i < freqYnUlZk.Count; i++)
                {
                    probYnUlZk = (double)freqYnUlZk[i] / (y.Count + u.Count + z.Count);
                    probYn1Zk = (double)freqYn1Zk[i] / (y.Count + z.Count);
                    probYnZk = (double)freqYnZk[i] / (y.Count + z.Count);
                    probYn1UlZk = (double)freqYn1UlZk[i] / (y.Count + u.Count + z.Count);
                    //probYnUl = (double)freqYnUl[i] / (y.Count + u.Count);
                    //probYn1Ul = (double)freqYn1Ul[i] / (y.Count + u.Count);
                    probYn1 = (double)freqYn1[i] / y.Count;
                    probYn = (double)freqYn[i] / y.Count;
                    if (Equals(probYn1UlZk, 0.0) || Equals(probYnUlZk, 0.0) || Equals(probYn1Zk, 0.0) || Equals(probYnZk, 0.0))
                        continue;
                    Hyn1 -= probYn1 * Math.Log(probYn1 / probYn, 2);
                    Hyn1ulzk -= probYn1UlZk * Math.Log(probYn1UlZk / probYnUlZk, 2);
                    DTE[lamda - 1] += probYn1UlZk * Math.Log(probYn1UlZk / probYnUlZk / (probYn1Zk / probYnZk), 2);
                }
                if (coef.Coeficiente < DTE[lamda - 1])
                {
                    coef.Coeficiente = DTE[lamda - 1];
                    coef.Lambda = lamda;
                }

                Y.RemoveAt(0);
                U.RemoveAt(U.Count - 1);
            }
            //coef.Coeficiente = coef.Coeficiente/Hy;
            return coef;
        }

        public static bool IsStacionary(RVector X, int m)
        {
            //Mean Test
            var x = new RMatrix();
            for (int i = 0; i < m; i++)
            {
                var temp = new List<double>();
                for (int j = 0; j < X.Length/m; j++)
                {
                    temp.Add(X[j]);
                }
                x.InsertColVector(new RVector(temp.ToArray()), i);
            }

            var mean = 0.0;
            for (int i = 0; i < m; i++)
            {
                mean += Summary.Mean(x.GetColVector(i));
            }
            mean = mean/m;

            var sigma = 0.0;
            for (int i = 0; i < m; i++)
            {
                var mj = Summary.Mean(x.GetColVector(i));
                sigma += Math.Pow(mj - mean, 2);
            }
            sigma = Math.Sqrt(sigma / (m*(m - 1)));

            for (int i = 0; i < m; i++)
            {
                var mj = Summary.Mean(x.GetColVector(i));
                if (mj < (mean - 6* sigma) && mj > (mean + 6*sigma))
                {
                    return false;
                }
            }
            return true;
        }

        public static RMatrix Histogram(List<double> y)
        {
            var codeBook = new List<double>();
            var min = y.Min();
            var max = y.Max();
            var h = (max - min)/10;
            foreach (var t in y.Where(t => codeBook.Contains(t)))
            {
                codeBook.Add(t);
            }
            codeBook.Sort();
            var frequency = codeBook.Select(t => y.Count(n => Equals(t, n))).Select(count => (double) count).ToList();

            var r = new RMatrix();
            r.InsertRowVector(new RVector(frequency.ToArray()), 0);
            r.InsertRowVector(new RVector(codeBook.ToArray()), 1);

            return r;
        }

        public static Correlacion NTE(IReadOnlyList<double> u, IReadOnlyList<double> y)
        {
            //var f = LloydAlgorithm(u, y);
            if (u == null) throw new ArgumentNullException(nameof(u));
            if (y == null) throw new ArgumentNullException(nameof(y));

            var coef = new Correlacion { Coeficiente = 0 };

            //var yHistogram = Summary.Histogram(new RVector(y.ToArray()));
            //var uHistogram = Summary.Histogram(new RVector(u.ToArray()));

            //var freqY = new List<double>();
            //var codeBookY = new List<double>();
            //for (int i = 0; i < yHistogram.NCols; i++)
            //{
            //    freqY.Add(yHistogram[0, i]);
            //    codeBookY.Add(yHistogram[1,i]);
            //}
            //var freqU = new List<double>();
            //var codeBookU = new List<double>();
            //for (int i = 0; i < uHistogram.NCols; i++)
            //{
            //    freqU.Add(uHistogram[0,i]);
            //    codeBookU.Add(uHistogram[1,i]);
            //}

            //var codeBook = new List<double>();
            //var min = y.Min();
            //var max = y.Max();
            //var h = (max - min) / 5;
            //codeBook.Add(min);
            //for (int i = 0; i < 5; i++)
            //{
            //    codeBook.Add(Math.Round(codeBook[i] + h, 2));
            //}
            //var Y = Cuntify(y, codeBook).ToList();
            //var U = Cuntify(u, codeBook).ToList();

            var len = 4;
            var Y = Cuntify(y, len).ToList();
            var U = Cuntify(u, len).ToList();

            //var n = UnivarModelOrder(y);
            //var n1 = n + 1;
            //var X = new RMatrix();
            //X.InsertColVector(Y, 0);
            //X.InsertColVector(U, 1);
            //var l = MultiModelOrder(0, X);
            //var l = UnivarModelOrder(U);
            var n = 1;
            var n1 = n + 1;
            var l = n;
            var NTE = new double[30];

            var codeBook = new List<double> { Y[0] };

            foreach (var t in Y.ToList().Where(t => !codeBook.Contains(t)))
                codeBook.Add(t);

            foreach (var t in U.ToList().Where(t => !codeBook.Contains(t)))
                codeBook.Add(t);

            codeBook.Sort();

            var freqY = new List<int>(codeBook.Count);
            var freqU = new List<int>(codeBook.Count);

            freqY.AddRange(codeBook.Select(t1 => Y.Count(t => Equals(t1, t))).Select(count => count));

            freqU.AddRange(codeBook.Select(t1 => U.Count(t => Equals(t1, t))).Select(count => count));

            var Hy = freqY.Where(t => !Equals(t, 0)).Aggregate(0.0, (current, t) => current - (double)t / Y.Count * Math.Log((double)t / Y.Count, 2));
            var Hu = freqU.Where(t => !Equals(t, 0)).Aggregate(0.0, (current, t) => current - (double)t / U.Count * Math.Log((double)t / U.Count, 2));

            for (var lamda = 0; lamda < 30; lamda++)
            {
                

                var freqYn1 = new List<int>();
                var freqYn1Ul = new List<int>();
                var freqYnUl = new List<int>();
                var freqYn = new List<int>();
                var freqUl = new List<int>();
                /*********************************************/
                foreach (var ty in codeBook)
                {
                    //foreach (var tu in codeBook)
                    //{
                        
                    //}
                    var cdY = new double[n];
                    var cdU = new double[l];
                    for (var i = 0; i < cdY.Length; i++)
                        cdY[i] = ty;
                    for (var i = 0; i < cdU.Length; i++)
                        cdU[i] = ty;

                    var count = 0;

                    for (int i = 0, m = 0; (i < Y.Count - n && m < U.Count - l);)
                    {
                        var tempY = new double[n];
                        for (var j = 0; j < n; j++)
                            tempY[j] = Y[i + j];

                        var tempU = new double[l];
                        for (var j = 0; j < l; j++)
                            tempU[j] = U[i + j];

                        if (Equals(new RVector(tempY), new RVector(cdY))
                            && Equals(new RVector(tempU), new RVector(cdU)))
                        {
                            count++;
                        }
                        i++;
                        m++;
                    }

                    freqYnUl.Add(count);
                }
                /*********************************************/
                foreach (var ty in codeBook)
                {
                    //foreach (var tu in codeBook)
                    //{
                        
                    //}
                    var cdYn1 = new double[n1];
                    var cdU = new double[l];
                    for (var i = 0; i < cdYn1.Length; i++)
                        cdYn1[i] = ty;
                    for (var i = 0; i < cdU.Length; i++)
                        cdU[i] = ty;

                    var count = 0;

                    for (int i = 0, m = 0; (i < Y.Count - n1 && m < U.Count - l);)
                    {
                        var tempYn1 = new double[n1];
                        for (var j = 0; j < n1; j++)
                            tempYn1[j] = Y[i + j];

                        var tempU = new double[l];
                        for (var j = 0; j < l; j++)
                            tempU[j] = U[i + j];

                        if (Equals(new RVector(tempYn1), new RVector(cdYn1))
                            && Equals(new RVector(tempU), new RVector(cdU)))
                        {
                            count++;
                        }
                        i++;
                        m++;
                    }

                    freqYn1Ul.Add(count);
                }
                /*********************************************/
                foreach (var ty in codeBook)
                {
                    //foreach (var tu in codeBook)
                    //{
                        
                    //}
                    var cdYn1 = new double[n1];
                    for (var i = 0; i < cdYn1.Length; i++)
                    {
                        cdYn1[i] = ty;
                    }
                    var count = 0;

                    for (int i = 0; i < Y.Count - n1;)
                    {
                        var temp = new double[n1];
                        for (var j = 0; j < n1; j++)
                        {
                            temp[j] = Y[i + j];
                        }
                        if (Equals(new RVector(temp), new RVector(cdYn1)))
                        {
                            count++;
                        }
                        i++;
                    }
                    freqYn1.Add(count);
                }
                /*********************************************/
                foreach (var ty in codeBook)
                {
                    //foreach (var tu in codeBook)
                    //{
                        
                    //}
                    var cdYn = new double[n];
                    for (var i = 0; i < cdYn.Length; i++)
                    {
                        cdYn[i] = ty;
                    }
                    var count = 0;

                    for (int i = 0; i < Y.Count - n;)
                    {
                        var temp = new double[n];
                        for (var j = 0; j < n; j++)
                        {
                            temp[j] = Y[i + j];
                        }
                        if (Equals(new RVector(temp), new RVector(cdYn)))
                        {
                            count++;
                        }
                        i++;
                    }
                    freqYn.Add(count);
                }
                /*********************************************/
                //foreach (var t in codeBook)
                //{
                //    var cdUn = new double[l];
                //    for (var i = 0; i < cdUn.Length; i++)
                //    {
                //        cdUn[i] = t;
                //    }
                //    var count = 0;

                //    for (int i = 1; i < U.Length / l; i++)
                //    {
                //        var temp = new double[l];
                //        for (var j = 0; j < l; j++)
                //        {
                //            temp[j] = Y[i + j];
                //        }
                //        if (Equals(new RVector(temp), new RVector(cdUn)))
                //        {
                //            count++;
                //        }
                //        i += l;
                //    }
                //    freqUl.Add(count);
                //}

                //var Hyn = freqYn.Where(t => !Equals(t, 0)).Aggregate(0.0, (current, t) => current - (double)t / Y.ToList().Count * Math.Log((double)t / Y.ToList().Count, 2));
                //var Hyn1 = freqYn1.Where(t => !Equals(t, 0)).Aggregate(0.0, (current, t) => current - (double)t / Y.ToList().Count * Math.Log((double)t / Y.ToList().Count, 2));
                
                var Hyn1 = 0.0;
                /*********************************************/
                var probYn1Ul = 0.0;
                var probYn1 = 0.0;
                var probYn = 0.0;
                var probYnUl = 0.0;
                /*********************************************/
                for (int i = 0; i < freqYn1.Count; i++)
                {
                    probYnUl = (double)freqYnUl[i] / (y.Count + u.Count);
                    probYn1 = (double)freqYn1[i] / y.Count;
                    probYn = (double)freqYn[i] / y.Count;
                    probYn1Ul = (double)freqYn1Ul[i] / (y.Count + u.Count);
                    if (Equals(probYn1Ul, 0.0) || Equals(probYnUl, 0.0) || Equals(probYn1, 0.0) || Equals(probYn, 0.0))
                        continue;
                    Hyn1 -= probYn1 * Math.Log(probYn1 / probYn, 2);
                    NTE[lamda] += probYn1Ul * Math.Log(probYn1Ul / probYnUl / (probYn1 / probYn), 2);
                }
                if (Equals(Hy, 0.0) /*|| Equals(Hyn1, 0.0) */|| Equals(NTE[lamda], 0.0))
                    continue;

                if (coef.Coeficiente < NTE[lamda])
                {
                    coef.Coeficiente = NTE[lamda];
                    coef.Lambda = lamda;
                }

                Y.RemoveAt(0);
                U.RemoveAt(U.Count - 1);
            }
            coef.Coeficiente = coef.Coeficiente/Hy;
            return coef;
        }

        public static Correlacion TE2(IReadOnlyList<double> u, IReadOnlyList<double> y)
        {
            if (u == null) throw new ArgumentNullException(nameof(u));
            if (y == null) throw new ArgumentNullException(nameof(y));

            var coef = new Correlacion { Coeficiente = 0 };

            var lambda = (int)Math.Sqrt(y.Count) / 3;
            
            //var delta = 1/(4*Math.Pow(U.Length, 0.25)*Math.Sqrt(Math.PI*Math.Log(U.Length)));
            for (var k = 0; k < 1; k++)
            {
                var len = 4;
                var Y = Cuntify(y, len);
                var U = Cuntify(u, len);

                var tempTE = 0.0;
                var Hy = 0.0;
                var X = new RMatrix();
                X.InsertColVector(Y, 0);
                X.InsertColVector(U, 1);
                Hy = Causality.Hy(Y);
                var modelOrderY = UnivarModelOrder(Y.ToList());
                var modelOrderYU = MultiModelOrder(0, X);

                /**********************************************/
                var varParamYnU = VarTEURM(X, modelOrderYU);
                var covYU = new RMatrix(varParamYnU.Model.NCols, varParamYnU.Model.NCols);
                for (int j = 0; j < varParamYnU.Model.NCols; j++)
                {
                    for (int l = 0; l < varParamYnU.Model.NCols; l++)
                    {
                        var vec = new RMatrix();
                        vec.InsertColVector(varParamYnU.Model.GetColVector(j), 0);
                        vec.InsertColVector(varParamYnU.Model.GetColVector(l), 1);
                        var temp = Summary.Covariance(vec);
                        if (j < l)
                            covYU[j, l] = temp[0, 1];
                        else
                            covYU[j, l] = temp[1, 0];
                    }
                }
                /**********************************************/
                var varParamYn = VarTERM(Y.ToList(), modelOrderY);
                var covY = new RMatrix(varParamYn.Model.NCols, varParamYn.Model.NCols);
                for (int j = 0; j < varParamYn.Model.NCols; j++)
                {
                    for (int l = 0; l < varParamYn.Model.NCols; l++)
                    {
                        var vec = new RMatrix();
                        vec.InsertColVector(varParamYn.Model.GetColVector(j), 0);
                        vec.InsertColVector(varParamYn.Model.GetColVector(l), 1);
                        var temp = Summary.Covariance(vec);
                        if (j < l)
                            covY[j, l] = temp[0, 1];
                        else
                            covY[j, l] = temp[1, 0];
                    }
                }
                /*********************************************/
                var varParam = JointPDF(X, modelOrderYU);
                var covYnYU = new RMatrix(varParam.Model.NCols, varParam.Model.NCols);
                for (int j = 0; j < varParam.Model.NCols; j++)
                {
                    for (int l = 0; l < varParam.Model.NCols; l++)
                    {
                        var vec = new RMatrix();
                        vec.InsertColVector(varParam.Model.GetColVector(j), 0);
                        vec.InsertColVector(varParam.Model.GetColVector(l), 1);
                        var temp = Summary.Covariance(vec);
                        if (j < l)
                            covYnYU[j, l] = temp[0, 1];
                        else
                            covYnYU[j, l] = temp[1, 0];
                    }
                }
                /**********************************************/

                var probJoint = 0.0;
                var probCondYnYU = 0.0;
                var probCondYnY = 0.0;
                for (int i = 0; i < varParam.Model.NRows; i++)
                {
                    var v = new RMatrix();
                    v.InsertColVector(varParam.Model.GetRowVector(i), 0);
                    var det = RMatrix.QRDeterminant(RMatrix.Inverse(covYnYU));
                    var temp = v.GetTranspose() * RMatrix.Inverse(covYnYU) * v;
                    probJoint = Math.Pow(2 * Math.PI, -1.0 * modelOrderYU / 2)
                        * Math.Pow(Math.Pow(varParam.Sigma, -2), -1.0 * modelOrderYU / 2)
                        * Math.Sqrt(det)
                        * Math.Exp(-0.5 * temp[0, 0] / Math.Pow(varParam.Sigma, 2))
                        ;
                    if (Equals(probJoint, 0.0))
                        continue;

                    probCondYnYU = Math.Exp(-0.5 * Math.Pow(varParamYnU.Epsilon[i] / varParamYnU.Sigma, 2))
                        / Math.Sqrt(2 * Math.PI * Math.Pow(varParamYnU.Sigma, 2));

                    probCondYnY = Math.Exp(-0.5 * Math.Pow(varParamYn.Epsilon[i] / varParamYn.Sigma, 2))
                        / Math.Sqrt(2 * Math.PI * Math.Pow(varParamYn.Sigma, 2));

                    tempTE += probJoint * Math.Log(probCondYnYU / probCondYnY, 2);
                }

                //for (var i = 2; i < 4; i++)
                //{
                //    var probYnYU = JointPDF(X, i + 1);
                //    var probYnU = VarTEURM(X, i + 1);
                //    var probYn = VarTERM(y, i);

                //    tempTE += probYnYU * Math.Log(probYnU / probYn, 2);
                //}
                //if (coef.Coeficiente < tempTE)
                //{
                //    coef.Coeficiente = tempTE;
                //    coef.Lambda = k;
                //}
                coef.Coeficiente = tempTE;
                coef.Lambda = k;
            }

            return coef;
            
        }

        private static Correlacion TE(IReadOnlyList<double> u, IReadOnlyList<double> y)
        {
            if (u == null) throw new ArgumentNullException(nameof(u));
            if (y == null) throw new ArgumentNullException(nameof(y));

            var coef = new Correlacion {Coeficiente = 0};

            
            

            var lambda = (int)Math.Sqrt(y.Count)/3;
            for (var j = 1; j < lambda; j++)
            {
                //var histogram = new RMatrix(Summary.Histogram(new RVector(y.ToArray()), (int)Math.Sqrt(y.Count)));
                //var Y = new RVector(histogram.NCols - j);
                //for (var i = 0; i < histogram.NCols - j; i++)
                //{
                //    Y[i] = histogram[0, i];
                //}

                //histogram = new RMatrix(Summary.Histogram(new RVector(u.ToArray()), (int)Math.Sqrt(u.Count)));
                //var U = new RVector(histogram.NCols - j);
                //for (var i = 0; i < histogram.NCols - j; i++)
                //{
                //    U[i] = histogram[0, i];
                //}

                //histogram = new RMatrix(Summary.Histogram(new RVector(y.ToArray()), (int)Math.Sqrt(y.Count)));
                //var Yn = new RVector(histogram.NCols - j);
                //for (var i = 0; i < histogram.NCols - j; i++)
                //{
                //    Yn[i] = histogram[0, i + j];
                //}
                var len = 4;
                var Y = Cuntify(y, len);
                var U = Cuntify(u, len);
                var Yn = Cuntify(y, len);
                
                //var tempTE = 0.0;
                //var Hy = 0.0;

                //var X = new RMatrix();
                //X.InsertColVector(Y, 0);
                //X.InsertColVector(U, 1);

                //Hy = Causality.Hy(Y);

                //var probYnYU = JointPDF(X, 2);
                //var probYnU = VarTEURM(X, 2);
                //var probYn = VarTERM(Y.ToList(), 2);

                var probYnYU = 0.0;
                var probYnY = 0.0;
                var probYU = 0.0;
                var probY = 0.0;
                var tempTE = 0.0;
                var Hy = 0.0;

                var YU = new RMatrix();
                var YnY = new RMatrix();
                var YnYU = new RMatrix();

                YU.InsertColVector(Y, 0);
                YU.InsertColVector(U, 1);

                YnY.InsertColVector(Yn, 0);
                YnY.InsertColVector(Y, 1);

                YnYU.InsertColVector(Yn, 0);
                YnYU.InsertColVector(Y, 1);
                YnYU.InsertColVector(U, 2);

                var covYU = Summary.Covariance(YU);
                var covYnY = Summary.Covariance(YnY);
                var covYnYU = Summary.Covariance(YnYU);

                if ((Math.Abs(RMatrix.QRDeterminant(covYU)) > 1e-15)
                    || (Math.Abs(RMatrix.QRDeterminant(covYnY)) > 1e-15)
                    || (Math.Abs(RMatrix.QRDeterminant(covYnYU)) > 1e-15))
                {
                    var uMedia = Summary.Mean(U);
                    var yMedia = Summary.Mean(Y);
                    var ynMedia = Summary.Mean(Yn);

                    var uSigma = Summary.StandardDeviation(U);
                    var ySigma = Summary.StandardDeviation(Y);
                    var ynSigma = Summary.StandardDeviation(Yn);

                    var sY = (int)Math.Round(1.06 * Math.Pow(Y.Length, -1.0 / 7.0) * ySigma);
                    var sU = (int)Math.Round(1.06 * Math.Pow(Y.Length, -1.0 / 6.0) * uSigma);
                    var sYn = (int)Math.Round(1.06 * Math.Pow(Y.Length, -0.2) * ynSigma);

                    var VecYU = new RMatrix(new double[2, 2]);
                    var VecYnY = new RMatrix(new double[2, 2]);
                    var VecYnYU = new RMatrix(new double[3, 3]);
                    var temp = new RMatrix(1, 1);

                    for (var i = 0; i < U.Length; i++)
                    {
                        double[] vecYnYU = { Yn[i] - ynMedia, Y[i] - yMedia, U[i] - uMedia };
                        VecYnYU = TransformCol(new RVector(vecYnYU));
                        double[] vecYnY = { Yn[i] - ynMedia, Y[i] - yMedia };
                        VecYnY = TransformCol(new RVector(vecYnY));
                        double[] vecYU = { Y[i] - yMedia, U[i] - uMedia };
                        VecYU = TransformCol(new RVector(vecYU));

                        temp = (VecYnYU.GetTranspose() * RMatrix.Inverse(covYnYU)) * VecYnYU;
                        probYnYU += (Math.Exp(-0.5 * temp[0, 0])
                            / Math.Sqrt(Math.Pow(2 * Math.PI, 3) * RMatrix.QRDeterminant(covYnYU)))
                            ;
                        temp = (VecYU.GetTranspose() * RMatrix.Inverse(covYU)) * VecYU;
                        probYU += (Math.Exp(-0.5 * temp[0, 0])
                            / Math.Sqrt(Math.Pow(2 * Math.PI, 2) * RMatrix.QRDeterminant(covYU)))
                            ;
                        temp = (VecYnY.GetTranspose() * RMatrix.Inverse(covYnY)) * VecYnY;
                        probYnY += (Math.Exp(-0.5 * temp[0, 0])
                            / Math.Sqrt(Math.Pow(2 * Math.PI, 2) * RMatrix.QRDeterminant(covYnY)))
                            ;

                        probY += Gauss(sY, Y[i]);

                        Hy -= probY * Math.Log(probY, 2);
                        
                    }

                    probY = probY / (Y.Length * sY);
                    probYU = probYU / (Y.Length * sY * sU);
                    probYnY = probYnY / (Y.Length * sY * sYn);
                    probYnYU = probYnYU / (Y.Length * sY * sYn * sU);
                    
                }


                tempTE = probYnYU*Math.Log(probYnYU/probYU
                                            /(probYnY/probY), 2);
                                            //-probYnY * Math.Log(probYnY / probY / probY, 2);


                if (Math.Abs(coef.Coeficiente) < Math.Abs(tempTE))
                {
                    coef.Coeficiente = tempTE;
                    coef.Lambda = j;
                }

            }

            return coef;
        }

        
        public static double Hy(RVector Y)
        {
            //var Y = new RVector(y.ToArray());
            var yMedia = Summary.Mean(Y);
            var ySigma = Summary.StandardDeviation(Y);
            
            var Hy = 0.0;
            for (int i = 0; i < Y.Length; i++)
            {
                var s = (int)Math.Round(1.06 * Math.Pow(Y.Length, -0.2) * ySigma);
                var probY = Gauss(s, Y[i]);

                if (!double.IsNaN(probY) && Math.Abs(probY) > 1e-15)
                    Hy -= probY*Math.Log(probY, 2);
            }
            return Hy;
        }

        public static double Gauss(double Sigma, double rss)
        {
            var a = 0.0;
            for (var i = -Sigma; i <= Sigma; i++)
            {
                a += Math.Exp(-1.0 * Math.Pow(i, 2) / Math.Pow(Sigma, 2)) / (Sigma * Math.Sqrt(Math.PI));
            }

            return Math.Exp(-0.5 * Math.Pow(rss, 2)/Sigma) / (Math.Sqrt(2 * Math.PI) * Math.Sqrt(Sigma));
        }

        public static Correlacion CCF(RVector X, RVector Y, int N, int cantidad)
        {
            if (X == null) throw new ArgumentNullException(nameof(X));
            if (Y == null) throw new ArgumentNullException(nameof(Y));

            var k1 = X.Length/4;
            var CoefCorr = new Correlacion();

            //var xMedia = Summary.Mean(X);
            //var yMedia = Summary.Mean(Y);

            
            var xDes = 0.0;
            var yDes = 0.0;
            for (var i = 0; i < k1; i++)
            {
                var xMedia = 0.0;
                var yMedia = 0.0;
                var covXY = 0.0;
                var tempCoef = 0.0;
                for (var j = 0; j < N - i; j++)
                {
                    xMedia += X[j];
                    yMedia += Y[j];
                }
                xMedia = xMedia/(N - i);
                yMedia = yMedia/(N - i);
                for (var j = 0; j < N - i; j++)
                {
                    covXY += (X[j] - xMedia) * (Y[j + i] - yMedia);
                    xDes += Math.Pow(X[j] - xMedia, 2);
                }
                xDes = xDes/(N - i);
                for (var j = i; j < N; j++)
                {
                    yDes += Math.Pow(Y[j] - yMedia, 2);
                }
                yDes = yDes/(N - i);

                xDes = Math.Sqrt(xDes);
                yDes = Math.Sqrt(yDes);

                //xDes = Summary.StandardDeviation(X);
                //yDes = Summary.StandardDeviation(Y);

                tempCoef = covXY / (xDes * yDes * X.Length);
                
                if (tempCoef > CoefCorr.Coeficiente)
                {

                    CoefCorr.Coeficiente = tempCoef;
                    CoefCorr.Lambda = i;
                }
            }
           
            return CoefCorr;
        }

        public static void ObtainCausality()
        {
            try
            {
                using (var sr = new StreamReader(AppDomain.CurrentDomain.BaseDirectory + "/CSTHma.txt"))
                {
                    var X = new RMatrix();
                    var lines = sr.ReadLine();
                    var names = new List<string>();
                    if (lines != null)
                    {
                        var temp = lines.Split('\t');
                        names = temp.ToList();
                        var valores = new double[temp.Length];

                        //Se capturan los datos historicos del proceso
                        while ((lines = sr.ReadLine()) != null)
                        {
                            temp = lines.Split('\t');
                            for (var i = 0; i < temp.Length; i++)
                            {
                                valores[i] = double.Parse(temp[i]);
                            }

                            X.InsertRowVector(new RVector(valores), X.NRows);
                        }
                    }
                    var stacionary = true;
                    for (int i = 0; i < X.NCols; i++)
                        if (!IsStacionary(X.GetColVector(i), 10))
                            stacionary = false;
                    
                    //Corss Correlation Funtion (CCF)
                    //Matriz con los Coeficientes de Correlacion Lineal de Pearson
                    var CoefCCF = new double[X.NCols, X.NCols];
                    //Matriz de Retardo
                    var RetardoCCF = new int[X.NCols, X.NCols];

                    var CoefCorrCCF = new Correlacion(0, 0);
                    for (int i = 0; i < X.NCols; i++)
                    {
                        for (int j = 0; j < X.NCols; j++)
                        {
                            if (i != j)
                            {
                                CoefCorrCCF = CCF(X.GetColVector(i), X.GetColVector(j), X.NRows, X.NCols);
                                RetardoCCF[i, j] = CoefCorrCCF.Lambda;
                                CoefCCF[i, j] = CoefCorrCCF.Coeficiente;

                            }
                            else
                                CoefCCF[i, j] = 0;
                        }
                    }

                    var Ro = new double[X.NCols, X.NCols];

                    for (int i = 0; i < X.NCols; i++)
                    {
                        for (int j = 0; j < X.NCols; j++)
                        {
                            if (i != j)
                                Ro[i, j] = CoefCCF[i, j] < CoefCCF[j, i] ? CoefCCF[j, i] : CoefCCF[i, j];
                        }
                    }

                    var N = 10000;
                    var Cccf = new double[X.NCols, X.NCols];

                    for (int i = 0; i < X.NCols; i++)
                    {
                        for (int j = 0; j < X.NCols; j++)
                        {
                            //Statistical test (3-Sigma)
                            if (Ro[i, j] >
                                1.85 * Math.Pow(X.NRows, -0.41) + 2.37 * Math.Pow(X.NRows, -0.53))
                            {
                                if (i != j && Ro[i, j] * Math.Sqrt((N - 2) / (1 - Math.Pow(Ro[i, j], 2))) > 1.697)
                                {
                                    //Directinaly test
                                    Cccf[i, j] = (CoefCCF[i, j] - CoefCCF[j, i]) / (CoefCCF[i, j] + CoefCCF[j, i]);
                                }
                            }
                            else
                                Cccf[i, j] = 0;
                            if (Cccf[i, j] < 0)
                            {
                                Cccf[i, j] = 0;
                            }
                        }
                    }
                    //var names = new List<string> { "CWFlow", "Level", "OutFlow", "TempOF", "Steam" };
                    Visualizar(new RMatrix(Cccf), names, "CCF");
                    
                    var CausalMatrix = new double[X.NCols, X.NCols];
                    CM = CausalMatrix;
                    for (int i = 0; i < X.NCols; i++)
                        for (int j = 0; j < X.NCols; j++)
                        {
                            //if (Equals(Cccf[i, j], 0.0)) continue;
                            //CM[i, j] += 1;

                            CM[i, j] += Cccf[i, j];
                        }

                    var CCFRootCauseList = names.ToDictionary(t => t, t => 0.0);

                    for (int i = 0; i < names.Count; i++)
                    {
                        for (int j = 0; j < names.Count; j++)
                        {
                            if (!Equals(Cccf[i, j], 0.0))
                                CCFRootCauseList[names[i]] += Cccf[i, j];
                        }
                    }
                }

                using (var sr = new StreamReader(AppDomain.CurrentDomain.BaseDirectory + "/csthTE.txt"))
                {
                    var X = new RMatrix();
                    var lines = sr.ReadLine();
                    var names = new List<string>();
                    if (lines != null)
                    {
                        var temp = lines.Split('\t');
                        names = temp.ToList();
                        var valores = new double[temp.Length];

                        //Se capturan los datos historicos del proceso
                        while ((lines = sr.ReadLine()) != null)
                        {
                            temp = lines.Split('\t');
                            for (var i = 0; i < temp.Length; i++)
                            {
                                valores[i] = double.Parse(temp[i]);
                            }

                            X.InsertRowVector(new RVector(valores), X.NRows);
                        }
                    }
                    var stacionary = true;
                    for (int i = 0; i < X.NCols; i++)
                        if (!IsStacionary(X.GetColVector(i), 10))
                            stacionary = false;


                    //var x = new RMatrix(X.NRows, X.NCols);
                    //for (int i = 0; i < X.NRows; i++)
                    //{
                    //    for (int j = 0; j < X.NCols; j++)
                    //    {
                    //        x[i, j] = X[X.NRows - i - 1, j];
                    //    }
                    //}

                    //var min = new RVector(X.NCols);
                    //var max = new RVector(X.NCols);
                    //for (int i = 0; i < X.NCols; i++)
                    //{
                    //    min[i] = 0;
                    //    max[i] = 1;
                    //}
                    //X.ZScoreNormalize();
                    //X.MinMaxNormalize(min, max);

                    //var TE = TransferEntropy(X);
                    //var DTE = DirectTransferEntropy(X, TE);
                    var min = new RVector(X.NCols);
                    var max = new RVector(X.NCols);
                    for (int i = 0; i < X.NCols; i++)
                    {
                        min[i] = 0;
                        max[i] = 1;
                    }
                    X.ZScoreNormalize();
                    X.MinMaxNormalize(min, max);
                    //var h = Summary.Histogram(X.GetColVector(0));
                    //var h1 = Summary.Histogram(X.GetColVector(1));
                    //var h2 = Summary.Histogram(X.GetColVector(2));
                    //var t1 = Task.Factory.StartNew(() => TransferEntropy(X));
                    var TE = TransferEntropy(X);
                    var ETE = new double[TE.NCols, TE.NCols];

                    for (int i = 0; i < TE.NCols; i++)
                        for (int j = 0; j < TE.NCols; j++)
                        {
                            ETE[i, j] = TE[i, j] - TE[j, i];
                            if (ETE[i, j] < 0.001)
                                ETE[i, j] = 0;
                        }

                    Visualizar(new RMatrix(ETE), names, "TE");

                    //var t2 = Task.Factory.StartNew(() => DirectTransferEntropy(X, TE));
                    var DTE = DirectTransferEntropy(X, new RMatrix(ETE));
                    Visualizar(new RMatrix(DTE), names, "DTE");

                    //var CausalMatrix = new double[X.NCols, X.NCols];
                    //CM = CausalMatrix;
                    for (int i = 0; i < X.NCols; i++)
                        for (int j = 0; j < X.NCols; j++)
                        {
                            //if (Equals(ETE[i, j], 0.0)) continue;
                            //CM[i, j] += 1;
                            //CM[i, j] += ETE[i, j] + DTE[i, j];
                            CM[i, j] += ETE[i, j];
                        }

                    //for (int i = 0; i < X.NCols; i++)
                    //    for (int j = 0; j < X.NCols; j++)
                    //    {
                    //        if (Equals(DTE[i, j], 0.0)) continue;
                    //        CM[i, j] += 1;
                    //    }
                    //Visualizar(new RMatrix(CausalMatrix), names, "final");
                    var TERootCauseList = names.ToDictionary(t => t, t => 0.0);
                    var DTERootCauseList = names.ToDictionary(t => t, t => 0.0);

                    for (int i = 0; i < names.Count; i++)
                    {
                        for (int j = 0; j < names.Count; j++)
                        {
                            if (!Equals(ETE[i, j], 0.0))
                                TERootCauseList[names[i]] += ETE[i, j];

                            if (!Equals(DTE[i, j], 0.0))
                                DTERootCauseList[names[i]] += DTE[i, j];
                        }
                    }
                }

                using (var sr = new StreamReader(AppDomain.CurrentDomain.BaseDirectory + "/csthGC.txt"))
                {
                    var X = new RMatrix();
                    var lines = sr.ReadLine();
                    var names = new List<string>();
                    if (lines != null)
                    {
                        var temp = lines.Split('\t');
                        names = temp.ToList();
                        var valores = new double[temp.Length];

                        //Se capturan los datos historicos del proceso
                        while ((lines = sr.ReadLine()) != null)
                        {
                            temp = lines.Split('\t');
                            for (var i = 0; i < temp.Length; i++)
                            {
                                valores[i] = double.Parse(temp[i]);
                            }

                            X.InsertRowVector(new RVector(valores), X.NRows);
                        }
                    }

                    var stacionary = true;
                    for (int i = 0; i < X.NCols; i++)
                        if (!IsStacionary(X.GetColVector(i), 10))
                            stacionary = false;
                    if (!stacionary)
                    {
                        Console.WriteLine("No es Estacionaria.");
                    }

                    //var t3 = Task.Factory.StartNew(() => GCausality(X));
                    //var t4 = Task.Factory.StartNew(() => GCausalityxyz(X));

                    //var inicio = DateTime.Now;
                    //var GC = GCausality(X);
                    //var GCxyz = GCausalityxyz(X);
                    //var fin = DateTime.Now;
                    //var tiempo = fin - inicio;

                    //var GC = new RMatrix();
                    ////var GCxy = new RMatrix();
                    //var task1 = Task.Factory.StartNew(() =>
                    //{
                    //    GC = GCausality(X);
                    //});

                    var GC = GCNew(X);
                    
                    //var task2 = Task.Factory.StartNew(() =>
                    //{
                    //    GCxy = GCausalityxy(X);
                    //});

                    //Task.WaitAll(task1);

                    //var EGCxy = new double[GCxy.NCols, GCxy.NCols];

                    //for (int i = 0; i < GCxy.NCols; i++)
                    //    for (int j = 0; j < GCxy.NCols; j++)
                    //    {
                    //        EGCxy[i, j] = GCxy[i, j] - GCxy[j, i];
                    //        if (EGCxy[i, j] < 0.0)
                    //            EGCxy[i, j] = 0;
                    //    }

                    Visualizar(GC, names, "GC");
                    //Visualizar(new RMatrix(EGCxy), names, "GCxyz");
                    
                    //var CausalMatrix = new double[X.NCols, X.NCols];
                    for (int i = 0; i < X.NCols; i++)
                        for (int j = 0; j < X.NCols; j++)
                        {
                            //CausalMatrix[i, j] = (GC[i, j] + GCxyz[i, j])/2;
                            //CM[i,j] += GC[i, j] + EGCxy[i, j];

                            //if (Equals(GC[i, j], 0.0)) continue;
                            //CM[i, j] += 1;

                            CM[i, j] += GC[i, j];
                            CM[i, j] = CM[i, j] / 3;
                            //CM[i, j] = Math.Round(CM[i, j], 4);
                        }
                    //for (int i = 0; i < X.NCols; i++)
                    //{
                    //    for (int j = 0; j < X.NCols; j++)
                    //    {
                    //        if (CM[i,j] < CM[j,i])
                    //        {
                    //            CM[i, j] = 0.0;
                    //        }
                    //    }
                    //}
                    Visualizar(new RMatrix(CM), names, "final");

                    var GCRootCauseList = names.ToDictionary(t => t, t => 0.0);
                    //var GCxyzRootCauseList = names.ToDictionary(t => t, t => 0.0);
                    var FinalRootCauseList = names.ToDictionary(t => t, t => 0.0);

                    for (int i = 0; i < names.Count; i++)
                    {
                        for (int j = 0; j < names.Count; j++)
                        {
                            if (!Equals(GC[i, j], 0.0))
                                GCRootCauseList[names[i]] += GC[i, j];

                            //if (!Equals(EGCxy[i, j], 0.0))
                            //    GCxyzRootCauseList[names[i]] += EGCxy[i, j];

                            if (!Equals(CM[i, j], 0.0))
                                FinalRootCauseList[names[i]] += CM[i, j];
                        }
                    }
                }
                
                


                Console.WriteLine();

                ///**************************************************************
                //Transfer Entropy (TE)
                //**************************************************************/
                ////Matrix con los Coeficientes de Correlacion Lineal de Pearson
                //var CoefTE = new double[X.NCols, X.NCols];
                ////Matrix de Retardo
                //var RetardoTE = new int[X.NCols, X.NCols];
                //var CoefCorrTE = new Correlacion(0, 0);
                //for (int i = 0; i < X.NCols; i++)
                //{
                //    for (int j = 0; j < X.NCols; j++)
                //    {
                //        if (i != j)
                //        {
                //            if (LeastSquaresTest(X.GetColVector(i), X.GetColVector(j)))
                //            {
                //                //CoefCorrTE = TransferEntropy(X.GetColVector(i), X.GetColVector(j), X.NRows);
                //                //CoefCorr = TransferEntropy(valores[i], valores[j], valores[i].Count); 
                //                CoefCorrTE = TE(X.GetColVector(i).ToList(), X.GetColVector(j).ToList());
                //                RetardoTE[i, j] = CoefCorrTE.lambda;
                //                CoefTE[i, j] = CoefCorrTE.coeficiente;
                //                //Pruba de Significacion Estadistica con un grado de confianza del 95%
                //                //1.96 / (Math.Sqrt(valores[i].Count - Retardo[i, j]
                //                if (Math.Abs(CoefTE[i, j]) < (1.85 * Math.Pow(X.NRows, -0.41) + (2.37 * Math.Pow(X.NRows, -0.53))))
                //                {
                //                    CoefTE[i, j] = 0;
                //                }
                //            }

                //        }
                //        else
                //            CoefTE[i, j] = 0;
                //    }
                //}
                Console.WriteLine(" hola");

                ////Se calcula el segundo test de significancia para hallar la direccionalidad
                //var Cccf = new double[valores.Count, valores.Count]; ;
                //for (int i = 0; i < valores.Count; i++)
                //{
                //    for (int j = 0; j < valores.Count; j++)
                //    {
                //        if (i != j)
                //        {
                //            Cccf[i, j] = (Math.Abs(CoefCCF[i, j]) - Math.Abs(CoefCCF[j, i])) / (Math.Abs(CoefCCF[i, j]) + Math.Abs(CoefCCF[j, i]));
                //        }
                //    }

                //}
                ////Se optiene el Peso de los Arcos del Grafo         
                //for (int i = 0; i < valores.Count; i++)
                //{
                //    for (int j = 0; j < valores.Count; j++)
                //    {
                //        if (i != j)
                //        {
                //            if (Cccf[i, j] >= (1.85 * Math.Pow(valores[i].Count, -0.41) + (2.37 * Math.Pow(valores[i].Count, -0.53))))
                //            {
                //                causalMatrix[i, j] = 1;
                //            }
                //            else if (Cccf[i, j] <= -(1.85 * Math.Pow(valores[i].Count, -0.41) + (2.37 * Math.Pow(valores[i].Count, -0.53))))
                //            {
                //                causalMatrix[i, j] = -1;
                //            }
                //            else
                //                causalMatrix[i, j] = 0;
                //        }
                //    }
                //}
                
            }
            catch (Exception e)
            {
                Console.WriteLine(e.Message + e.StackTrace);
                
            }
        }
    }

    
    
}
