using System;
using System.Collections.Generic;
using System.Text;

namespace Optimum
{
    class Matrix
    {
        protected int rows, columns;
        protected double[,] data;
        public Matrix(int r, int c)
        {
            this.rows = r; this.columns = c;
            data = new double[rows, columns];
            for (int i = 0; i < rows; i++)
                for (int j = 0; j < columns; j++) data[i, j] = 0;
        }
        public Matrix(double[,] mm)
        {
            this.rows = mm.GetLength(0); this.columns = mm.GetLength(1);
            data = new double[rows, columns];
            for (int i = 0; i < rows; i++)
                for (int j = 0; j < columns; j++)
                    data[i, j] = mm[i, j];
        }
        public int GetCountRows()
        {
            return rows;
        }
        public int GetCountColumns()
        {
            return columns;
        }
        public double this[int i, int j]
        {
            get
            {
                if (i < 0 && j < 0 && i >= rows && j >= columns)
                {
                    // Console.WriteLine(" Индексы вышли за пределы матрицы ");
                    return Double.NaN;
                }
                else
                    return data[i, j];
            }
            set
            {
                if (i < 0 && j < 0 && i >= rows && j >= columns)
                {
                    //Console.WriteLine(" Индексы вышли за пределы матрицы ");
                }
                else
                    data[i, j] = value;
            }
        }
        public Vector GetRow(int r)
        {
            if (r >= 0 && r < rows)
            {
                Vector row = new Vector(columns);
                for (int j = 0; j < columns; j++) row[j] = data[r, j];
                return row;
            }
            return null;
        }
        public void SetRow(int r,Vector rr)
        {
            if (r >= 0 && r < rows)
            {
                if(columns==rr.Size)
                for (int j = 0; j < columns; j++)  data[r, j]= rr[j];
                
            }
        }
        public void SetColumn(int c, Vector cc)  //Набор столбцов
        {
            if (c >= 0 && c < columns)
            {
                if (rows == cc.Size)
                    for (int i = 0; i < rows; i++) data[i, c] = cc[i];
            }
        }
        public Vector GetColumn(int c)
        {
            if (c >= 0 && c < columns)
            {
                Vector column = new Vector(rows);
                for (int i = 0; i < rows; i++) column[i] = data[i, c];
                return column;
            }
            return null;
        }
        public static Vector operator *(Matrix a, Vector b)
        {
            if (a.columns != b.Size) return null;
            Vector r = new Vector(a.rows);
            for (int i = 0; i < a.rows; i++)
            {
                r[i] = a.GetRow(i) * b;
            }
            return r;
        }
        public void PrintMatrix()
        {
            for (int i = 0; i < rows; i++)
            {
                for (int j = 0; j < columns; j++)
                {
                    Console.Write($"{data[i, j]} \t");
                }
                Console.WriteLine();
            }
            Console.WriteLine("\n");
        }
        public Matrix Copy()
        {
            Matrix t = new Matrix(rows, columns);
            for (int i = 0; i < rows; i++)
            {
                for (int j = 0; j < columns; j++)
                {

                    t.data[i, j] = data[i, j];


                }

            }
            return t;
        }
        public override string ToString()
        {
            string s = "{";
            for (int i = 0; i < rows - 1; i++)
            {
                Vector v = this.GetRow(i);
                s += v.ToString() + ";";
            }
            Vector vc = this.GetRow(rows - 1);
            s += vc.ToString() + "}";
            return s;
        }
        public Matrix Transpr()//транспонирует матрицу
        {

            Matrix t = new Matrix(columns, rows);
            for (int i = 0; i < columns; i++)
            {
                for (int j = 0; j < rows; j++)
                {

                    t.data[i, j] = data[j, i];

                }
            }
            return t;
        }
        public static Matrix Umnch(Matrix a, double ch)//Умножаем матрицу на число
        {
            Matrix resMass = new Matrix(a.rows, a.columns);
            for (int i = 0; i < a.rows; i++)
            {
                for (int j = 0; j < a.columns; j++)
                {
                    resMass[i, j] = a[i, j] * ch;
                }
            }
            return resMass;
        }
        
        public static Matrix UmnMatrix(Matrix a, Matrix b)//Умножение матриц
        {
            if (a.GetCountColumns() != b.GetCountRows())
            {
                throw new Exception("Умножение не возможно");
            }

            Matrix rez = new Matrix(a.rows, b.columns);

            for (int i = 0; i < a.rows; i++)
            {
                for (int j = 0; j < b.columns; j++)
                {
                    rez[i, j]= 0;
                    for (int k = 0; k < a.columns; k++)
                    {
                        rez[i, j] = rez[i, j] + a[i, k] * b[k, j];
                    }
                }
            }
            return rez;
        }
        public static Matrix SumMatrix(Matrix a, Matrix b)//Сложение матриц
        {
            Matrix rez = new Matrix(a.rows, a.columns);

            for (int i = 0; i < a.rows; i++)
            {
                for (int j = 0; j < b.columns; j++)
                {
                    rez[i, j] += a[i, j] + b[i, j];
                }
            }
            return rez;
        }

        public static Matrix SubtractionMatrix(Matrix a, Matrix b)//Вычитание матриц
        {
            Matrix rez = new Matrix(a.rows, a.columns);

            for (int i = 0; i < a.rows; i++)
            {
                for (int j = 0; j < b.columns; j++)
                {
                    rez[i, j] += a[i, j] - b[i, j];
                }
            }
            return rez;
        }
        public static Vector SLU_DOWN(Matrix a, Vector b)
        {
            if (a.rows != b.Size) return null;
            for (int i = 0; i < a.rows; i++)
                for (int j = i + 1; j < a.columns; j++) if (a[i, j] != 0.0) return null;
            for (int i = 0; i < a.columns; i++) if (a[i, i] == 0.0) return null;
            Vector x = new Vector(b.Size);

            x[0] = b[0] / a[0, 0];
            x[1] = b[1] / (a[1, 1] + x[0] * a[1, 0]);
            x[2] = b[2] / (a[2, 2] + x[0] * a[2, 0] + x[1] * a[2, 1]);
            x[3] = b[3] / (a[3, 3] + x[0] * a[3, 0] + x[1] * a[3, 1] + x[2] * a[3, 2]);
            for (int i = 2; i < b.Size; i++)
            {
                double c = 0;
                for (int j = 0; j < i; j++) { c += x[j] * a[i, j]; }
                x[i] = (b[i] - c) / a[i, i];
            }
            return x;

        }
        public static Vector SLU_UP(Matrix a, Vector b)
        {
            if (a.rows != b.Size) return null;
            for (int i = a.rows - 1; i >= 0; i--)
                for (int j = i - 1; j >= 0; j--)
                {
                    if (a[i, j] != 0.0) return null;
                    if (a[i, i] == 0.0) return null;
                }
            Vector x = new Vector(b.Size);
            x[b.Size - 1] = b[b.Size - 1] / a[a.rows - 1, a.rows - 1];
            for (int i = a.rows - 2; i >= 0; i--)
            {
                double c = 0;
                for (int j = b.Size - 1; j > i; j--) { c += x[j] * a[i, j]; }
                x[i] = (b[i] - c) / (a[i, i]);
            }
            return x;

        }
        public static double OpredDet123(Matrix a)
        {
            if ((a.columns == a.rows) & (a.columns == 1)) { return a[0, 0]; }
            double s = 0;
            if ((a.columns == a.rows) & (a.columns == 3))
            {
                s = a[0, 0] * a[1, 1] * a[2, 2] + a[0, 1] * a[1, 2] * a[2, 0] + a[1, 0] * a[2, 1] * a[0, 2] - a[0, 2] * a[1, 1] * a[2, 0] - a[0, 1] * a[1, 0] * a[2, 2] - a[0, 0] * a[1, 2] * a[2, 1];
            }
            else if ((a.columns == a.rows) & (a.columns == 2))
            {
                s = a[0, 0] * a[1, 1] - a[0, 1] * a[1, 0];
            }
            else throw new Exception("Не 3 и не 2 ");
            return s;
        }
        public static Matrix Obratnay(Matrix a)
        {//ТОлько для 3
            if (a.columns != a.rows) return null;
            if (Matrix.OpredDet123(a) == 0) return null;
            var res = new Matrix(a.rows, a.columns);
            res = (1 / Math.Abs(Matrix.OpredDet123(a))) * a.Transpr();
            return res;

        }
        /*public static double Opredel(Matrix a)
         { 
             if (a.columns != a.rows) throw new Exception("Особая ");
             double s = 0;
             Matrix b = new Matrix(a.rows - 1, a.columns - 1);
             for (int k = 0; k < a.columns-1; k++)
             {
                 for (int i = 1; i < a.rows-1; i++)
                 {
                     for (int j = 0; j < a.columns-1; j++)
                     {
                         if (j != k)
                         {
                             b[i, j] = a[i, j];
                         }

                     }
                 }
                 if ((b.columns == 3) || (b.columns == 2)||(b.columns==1))
                 {
                     s += a[0, k] * OpredDet123(b);
                 }
                 else
                 {
                     s += a[0, k] * Opredel(b);
                 }


             }
             return s;
         }*/
        public static Vector Gauss(Matrix aa, Vector bb)
        {
            Matrix a = aa.Copy();
            Vector b = bb.Copy();
            Vector x;
            if (a.columns != b.Size) return null;
            if (a.columns != a.rows) return null;
            for (int j = 0; j < a.rows - 1; j++)
            {
                double maxc = 0;
                int jmax = j;
                for (int i = j; i < a.rows; i++)
                {
                    if (Math.Abs(a[i, j]) > maxc) { maxc = Math.Abs(a[i, j]); jmax = i; }
                }

                if (maxc == 0) return null;
                Vector rowj = a.GetRow(j);
                Vector rowmaxj = a.GetRow(jmax);
                for (int z = 0; z < a.columns; z++)//меняем местами столбики
                {
                    a[j, z] = rowmaxj[z];
                    a[jmax, z] = rowj[z];
                }
                double c = b[j]; b[j] = b[jmax]; b[jmax] = c;

                for (int i = j + 1; i < a.rows; i++)
                {
                    double t = a[i, j] / a[j, j];
                    for (int z = 0; z < a.columns; z++)
                    {
                        a[i, z] = a[i, z] - a[j, z] * t;
                    }
                    b[i] = b[i] - b[j] * t;
                }

            }
            x = SLU_UP(a, b);
            return x;
        }
        public static Vector MethodPoregonki(Vector v, Vector s, Vector n, Vector rav)//дебаг
        {
            int vSize = v.size;//верхняя диагональ
            int sSize = s.size;//средняя диагональ
            int nSize = n.size;//нижняя диагональ
            if (nSize + 1 != sSize || sSize != vSize + 1 || nSize != vSize)
            {
                return null;
            }

            Vector a = new Vector(sSize);//v
            Vector y = new Vector(sSize);//s
            Vector b = new Vector(sSize);//n

            y[0] = s[0];
            a[0] = -v[0] / y[0];
            b[0] = rav[0] / y[0];
            for (int i = 1; i < sSize-1 ; i++)
            {

                y[i] = s[i] + (n[i-1] * a[i - 1]);
                a[i] = -v[i] / y[i];
                b[i] = (rav[i] - n[i-1] * b[i - 1]) / y[i];
            }
            y[sSize - 1] = s[sSize - 1] + n[sSize - 2] * a[sSize - 2];

            b[sSize - 1] = (rav[sSize - 1] - n[sSize - 2] * b[sSize - 2]) / y[sSize - 1];//Возможно тут баг последний обнуляется 

            Vector X = new Vector(rav.Size);
            X[sSize - 1] = b[sSize - 1];
            for (int i = sSize - 2; i >= 0; i--)
            {
                X[i] = a[i] * X[i + 1] + b[i];
            }//либо тут 
            
            return X;
        }
       

        public static Vector MethodGrammaShmidta(Matrix A, Vector res)//353
        {
            var vec = new List<Vector>(A.rows);
            Matrix R = new Matrix(A.rows, A.rows);
            Matrix T = new Matrix(A.rows, A.rows);
            if (A.rows != A.columns) return null;
            
                for (int j = 0; j < A.columns; j++)
                {
                int p = 0;
                var z= new Vector(A.rows);
                for (int i = 0; i < A.rows; i++)
                    {
                    
                    double x = A[i, j];
                    z[p] = x;
                    p++;
                }
                    vec.Add(z);
                }
            
            //метод y
            for (int i = 0; i < R.rows; i++)
            {
                R[i, 0] = vec[0][i];

            }
            for (int i = 0; i < R.rows; i++)
            {
                
               var a = vec[i];
                for (int k = 1; k <= i; k++) {
                    var b = (vec[i] * vec[k - 1]) / (vec[k - 1] * vec[k - 1]) * vec[k - 1];//при последней итерации b должен быть 0.1714*век1
                     T[k-1,i] = (vec[i] * vec[k - 1]) / (vec[k - 1] * vec[k - 1]);
                    

                    a -= b;
                }
                /*(vec[i] * vec[i - 1]) /( vec[i - 1] * vec[i - 1])*vec[i-1]*/
                for (int j = 0; j < R.rows; j++)
                    {
                        R[j, i] = a[j];
                    }
                    
                
            }

            
            for (int i = 0; i < R.rows; i++)
            {
                for (int j = 0; j < R.columns; j++)
                {
                    if (i == j) T[i, j] = 1;
                    
                }
            }
            

            

        var R1 = R.Transpr();
            var b1 = R1 * res;
            var D = R1*R;

            for (int i = 0; i < D.rows; i++)
            {
                for (int j = 0; j < D.columns; j++)
                {
                    if (i == j) D[i, j] = 1 / D[i, j];
                    else D[i, j] = 0;
                    
                }

            }
           
            b1 = D * b1;
            Vector X = SLU_UP(T,b1) ;
            return X;
        }

        public static Vector MetodPoslPre(Matrix A, Vector B, double eps)
        {
            for (int i = 0; i < A.columns; i++)
            {
                if (A[i, i] == 0) return null;
            }
            if (A.rows != A.columns || B.Size != A.columns) return null;
            var a = new Matrix(A.rows, A.columns);
            var b = new Vector(B.Size);
            for (int k = 0; k < A.columns; k++)
            {
                for (int j = 0; j < A.rows; j++)
                {
                    if (k == j) a[k, j] = 0;
                    else a[k, j] = -A[k, j] / A[k, k];
                }
                b[k] = B[k] / A[k, k];
            }
            Vector x = b;
            double z = 0;
            Vector x1 = B;
            while (x1.NormaE() - z > eps)
            {
                x1 = b + a * x;
                z = x.NormaE();
                x = x1;

            }
            return x1;

        }





        public static Matrix operator *(Matrix a, Matrix b)
        {
            return Matrix.UmnMatrix(a, b);
        }
        

        public static Matrix operator *(Matrix a, double b)
        {
            return Matrix.Umnch(a, b);
        }

        public static Matrix operator *(double b, Matrix a)
        {
            return Matrix.Umnch(a, b);
        }

        public static Matrix operator +(Matrix a, Matrix b)
        {
            return Matrix.SumMatrix(a, b);
        }

        public static Matrix operator -(Matrix a, Matrix b)
        {
            return Matrix.SubtractionMatrix(a, b);
        }

    }
    }


