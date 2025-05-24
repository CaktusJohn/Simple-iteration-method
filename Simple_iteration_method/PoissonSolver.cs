using System;

namespace LAB_CHM_2023_3_1
{
    public class Solution
    {
        public class Result
        {
            public double[][] V1 { get; set; }
            public double[][] U { get; set; }
            public double[][] V2 { get; set; }
            public double[][] V2_2 { get; set; }
            public double[] X { get; set; }
            public double[] Y { get; set; }
            public int P { get; set; }
            public double EpsMax { get; set; }
            public double MaxPogr { get; set; }
            public double MaxF { get; set; }
            public double MaxR1 { get; set; }
            public double XMax { get; set; }
            public double YMax { get; set; }
            public int P2 { get; set; }
            public double EpsMax2 { get; set; }
            public double MaxF2 { get; set; }
            public double MaxR { get; set; }
        }

        public double U1(double x, double y) // U* Решение тестовой задачи
        {
            return Math.Exp(1 - Math.Pow(x, 2) - Math.Pow(y, 2));
        }

        public double F1(double x, double y) // Функция полученная через Лапласса
        {
            return -4 * Math.Exp(1 - Math.Pow(x, 2) - Math.Pow(y, 2)) * (x * x + y * y - 1);
        }

        public double F2(double x, double y) // F*
        {
            return Math.Abs(Math.Pow(x, 2) - Math.Pow(y, 2));
        }

        public double Mu1(double y) //Граничное условие 1
        {
            return -y * y + 1;
        }

        public double Mu2(double y) //Граничное условие 2
        {
            return (1 - y * y) * Math.Exp(y);
        }

        public double Mu3(double x) //Граничное условие 3
        {
            return 1 - x * x;
        }

        public double Mu4(double x) //Граничное условие 4
        {
            return 1 - x * x;
        }

        public Result SolveTestTask(int n, int m, int N_max, double Eps, double Tau)
        {
            double h = 2.0 / (double)n, k = 2.0 / (double)m; //Шаги по x и y
            double h2 = -1.0 / (h * h), k2 = -1.0 / (k * k);
            double A = -2 * (h2 + k2);
            double[][] f; //вектор правой части
            double[] x, y; //границы по х и по y
            double[][] R; // невязка
            int p = 0; //Текущее число итераций
            double MaxPogr = 0.0;
            double MaxF = 0.0;
            double maxR1 = 0.0;

            x = new double[n + 1];
            y = new double[m + 1];
            double[][] v1 = new double[n + 1][];
            f = new double[n + 1][];
            double[][] u = new double[n + 1][];
            R = new double[n + 1][];

            for (int i = 0; i <= n; i++)
            {
                v1[i] = new double[m + 1];
                f[i] = new double[m + 1];
                u[i] = new double[m + 1];
                R[i] = new double[m + 1];
            }

            for (int i = 0; i <= n; i++)  //Заполнение массива x
            {
                x[i] = -1 + i * h;
            }

            for (int j = 0; j <= m; j++)  //Заполнение массива y
            {
                y[j] = -1 + j * k;
            }

            for (int j = 0; j <= m; j++)            //Заполнение массивов f и u
            {
                for (int i = 0; i <= n; i++)
                {
                    f[i][j] = F1(x[i], y[j]);
                    u[i][j] = U1(x[i], y[j]);
                    if (Math.Abs(f[i][j]) > MaxF) MaxF = Math.Abs(f[i][j]);
                    R[i][j] = 0;
                }
            }

            for (int j = 0; j <= m; j++)  //Заполнение граничных условий в массив v1
            {
                v1[0][j] = U1(-1, y[j]);
                v1[n][j] = U1(1, y[j]);
            }

            for (int i = 0; i <= n; i++)  //Заполнение граничных условий в массив v1
            {
                v1[i][0] = U1(x[i], -1);
                v1[i][m] = U1(x[i], 1);
            }

            for (int j = 1; j < m; j++)    //Нулевое начальное приближение
            {
                for (int i = 1; i < n; i++)
                {
                    v1[i][j] = 0.0;
                }
            }

            // МЕТОД ПРОСТОЙ ИТЕРАЦИИ
            double temp, prev, currentEps;
            double Eps_max;

            while (true)
            {
                for (int j = 1; j < m; j++)
                {
                    for (int i = 1; i < n; i++)
                    {
                        R[i][j] = A * v1[i][j] + h2 * (v1[i - 1][j] + v1[i + 1][j]) + k2 * (v1[i][j - 1] + v1[i][j + 1]) - F1(x[i], y[j]);
                    }
                }

                Eps_max = 0.0;
                for (int j = 1; j < m; j++)
                {
                    for (int i = 1; i < n; i++)
                    {
                        prev = v1[i][j];
                        temp = prev - Tau * R[i][j];
                        currentEps = Math.Abs(prev - temp);
                        if (currentEps > Eps_max) { Eps_max = currentEps; };
                        v1[i][j] = temp;
                    }
                }

                p++;
                if ((Eps_max < Eps) || (p > N_max))
                    break;
            }

            temp = 0;
            for (int j = 1; j < m; j++)
            {
                for (int i = 1; i < n; i++)
                {
                    temp = A * v1[i][j] + h2 * (v1[i - 1][j] + v1[i + 1][j]) + k2 * (v1[i][j - 1] + v1[i][j + 1]) - F1(x[i], y[j]);
                    if (Math.Abs(temp) >= maxR1) maxR1 = Math.Abs(temp);
                }
            }

            double xMax = 0.0;
            double yMax = 0.0;
            double Pogr;
            for (int j = 0; j <= m; j++)
            {
                for (int i = 0; i <= n; i++)
                {
                    Pogr = Math.Abs(u[i][j] - v1[i][j]);
                    if (Pogr > MaxPogr)
                    {
                        MaxPogr = Pogr;
                        xMax = x[i];
                        yMax = y[j];
                    }
                }
            }

            return new Result
            {
                V1 = v1,
                U = u,
                X = x,
                Y = y,
                P = p,
                EpsMax = Eps_max,
                MaxPogr = MaxPogr,
                MaxF = MaxF,
                MaxR1 = maxR1,
                XMax = xMax,
                YMax = yMax
            };
        }

        public Result SolveMainTask(int n, int m, int N_max, double Eps, double Tau)
        {
            double h = 2.0 / n, k = 2.0 / m; //Шаги по x и y
            double h2 = -1.0 / (h * h), k2 = -1.0 / (k * k);
            double A = -2 * (h2 + k2);
            double[][] f;
            double[][] R1; //невязка
            double[][] R2; //невязка
            double[] x, y;

            int p = 0; //Текущее число итераций
            double MaxPogr = 0.0;
            double MaxF = 0.0, MaxF2 = 0.0;
            double maxR1 = 0.0;

            x = new double[n + 1];
            y = new double[m + 1];
            double[][] v2 = new double[n + 1][];
            f = new double[n + 1][];
            R1 = new double[n + 1][];

            for (int i = 0; i <= n; i++)
            {
                v2[i] = new double[m + 1];
                f[i] = new double[m + 1];
                R1[i] = new double[m + 1];
            }

            for (int i = 0; i <= n; i++)  //Заполнение массива x
            {
                x[i] = -1 + i * h;
            }

            for (int j = 0; j <= m; j++)  //Заполнение массива y
            {
                y[j] = -1 + j * k;
            }

            for (int j = 0; j <= m; j++)            //Заполнение массива f
            {
                for (int i = 0; i <= n; i++)
                {
                    f[i][j] = F2(x[i], y[j]);
                    if (Math.Abs(f[i][j]) > MaxF) MaxF = Math.Abs(f[i][j]);
                    R1[i][j] = 0;
                }
            }

            for (int j = 0; j <= m; j++)  //Заполнение граничных условий в массив v2
            {
                v2[0][j] = Mu1(y[j]);
                v2[n][j] = Mu2(y[j]);
            }

            for (int i = 0; i <= n; i++)  //Заполнение граничных условий в массив v2
            {
                v2[i][0] = Mu3(x[i]);
                v2[i][m] = Mu4(x[i]);
            }

            for (int j = 1; j < m; j++)    //Нулевое начальное приближение
            {
                for (int i = 1; i < n; i++)
                {
                    v2[i][j] = 0.0;
                }
            }

            // UpRelaxMethod
            double temp, prev, currentEps;
            double Eps_max;

            while (true)
            {
                for (int j = 1; j < m; j++)
                {
                    for (int i = 1; i < n; i++)
                    {
                        R1[i][j] = A * v2[i][j] + h2 * (v2[i - 1][j] + v2[i + 1][j]) + k2 * (v2[i][j - 1] + v2[i][j + 1]) - F2(x[i], y[j]);
                    }
                }

                Eps_max = 0.0;
                for (int j = 1; j < m; j++)
                {
                    for (int i = 1; i < n; i++)
                    {
                        prev = v2[i][j];
                        temp = prev - Tau * R1[i][j];
                        currentEps = Math.Abs(prev - temp);
                        if (currentEps > Eps_max) { Eps_max = currentEps; };
                        v2[i][j] = temp;
                    }
                }

                p++;
                if ((Eps_max < Eps) || (p > N_max))
                    break;
            }

            // nevyazka na vyhode
            temp = 0.0;
            for (int j = 1; j < m; j++)
            {
                for (int i = 1; i < n; i++)
                {
                    temp = A * v2[i][j] + h2 * (v2[i - 1][j] + v2[i + 1][j]) + k2 * (v2[i][j - 1] + v2[i][j + 1]) - F2(x[i], y[j]);
                    if (Math.Abs(temp) >= maxR1) maxR1 = Math.Abs(temp);
                }
            }

            // solution whis step / 2
            n = 2 * n;
            m = 2 * m;
            x = new double[n + 1];
            y = new double[m + 1];
            double[][] v2_2 = new double[n + 1][];
            R2 = new double[n + 1][];
            f = new double[n + 1][];

            h = 2.0 / n;
            k = 2.0 / m;
            h2 = -1.0 / (h * h);
            k2 = -1.0 / (k * k);
            A = -2 * (h2 + k2);

            int p2 = 0;
            double maxR = 0.0;

            for (int i = 0; i <= n; i++)
            {
                v2_2[i] = new double[m + 1];
                R2[i] = new double[m + 1];
                f[i] = new double[m + 1];
            }

            for (int i = 0; i <= n; i++)  //Заполнение массива x
            {
                x[i] = -1 + i * h;
            }

            for (int j = 0; j <= m; j++)  //Заполнение массива y
            {
                y[j] = -1 + j * k;
            }

            for (int j = 0; j <= m; j++)            //Заполнение массива f
            {
                for (int i = 0; i <= n; i++)
                {
                    f[i][j] = F2(x[i], y[j]);
                    if (Math.Abs(f[i][j]) > MaxF2) MaxF2 = Math.Abs(f[i][j]);
                }
            }

            for (int j = 0; j <= m; j++)  //Заполнение граничных условий в массив v
            {
                v2_2[0][j] = Mu1(y[j]);
                v2_2[n][j] = Mu2(y[j]);
            }

            for (int i = 0; i <= n; i++)  //Заполнение граничных условий в массив v
            {
                v2_2[i][0] = Mu3(x[i]);
                v2_2[i][m] = Mu4(x[i]);
            }

            for (int j = 1; j < m; j++)    //Нулевое начальное приближение
            {
                for (int i = 1; i < n; i++)
                {
                    v2_2[i][j] = 0.0;
                }
            }

            // UpRelaxMethod
            temp = 0.0;
            prev = 0.0;
            currentEps = 0.0;
            double Eps_max2;

            while (true)
            {
                for (int j = 1; j < m; j++)
                {
                    for (int i = 1; i < n; i++)
                    {
                        R2[i][j] = A * v2_2[i][j] + h2 * (v2_2[i - 1][j] + v2_2[i + 1][j]) + k2 * (v2_2[i][j - 1] + v2_2[i][j + 1]) - F2(x[i], y[j]);
                    }
                }

                Eps_max2 = 0.0;
                for (int j = 1; j < m; j++)
                {
                    for (int i = 1; i < n; i++)
                    {
                        prev = v2_2[i][j];
                        temp = prev - Tau * R2[i][j];
                        currentEps = Math.Abs(prev - temp);
                        if (currentEps > Eps_max2) { Eps_max2 = currentEps; };
                        v2_2[i][j] = temp;
                    }
                }

                p2++;
                if ((Eps_max2 < Eps) || (p2 > N_max))
                    break;
            }

            // nevyazka na vyhode
            temp = 0.0;
            for (int j = 1; j < m; j++)
            {
                for (int i = 1; i < n; i++)
                {
                    temp = A * v2_2[i][j] + h2 * (v2_2[i - 1][j] + v2_2[i + 1][j]) + k2 * (v2_2[i][j - 1] + v2_2[i][j + 1]) - F2(x[i], y[j]);
                    if (Math.Abs(temp) >= maxR) maxR = Math.Abs(temp);
                }
            }

            n = n / 2;
            m = m / 2;

            double xMax = 0.0;
            double yMax = 0.0;
            double Pogr;
            for (int j = 0; j <= m; j++)
            {
                for (int i = 0; i <= n; i++)
                {
                    Pogr = Math.Abs(v2[i][j] - v2_2[i * 2][j * 2]);
                    if (Pogr > MaxPogr)
                    {
                        MaxPogr = Pogr;
                        xMax = x[2 * i];
                        yMax = y[2 * j];
                    }
                }
            }

            return new Result
            {
                V2 = v2,
                V2_2 = v2_2,
                X = x,
                Y = y,
                P = p,
                EpsMax = Eps_max,
                MaxPogr = MaxPogr,
                MaxF = MaxF,
                MaxR1 = maxR1,
                P2 = p2,
                EpsMax2 = Eps_max2,
                MaxF2 = MaxF2,
                MaxR = maxR,
                XMax = xMax,
                YMax = yMax
            };
        }

        public double CalculateTau(int n, int m)
        {
            double h = 2.0 / (double)n, k = 2.0 / (double)m; //Шаги по x и y
            double h2 = -1.0 / (h * h), k2 = -1.0 / (k * k);
            double Max, Min;
            Min = -4 * h2 * Math.Pow(Math.Sin(Math.PI / (2.0 * n)), 2) - 4 * k2 * Math.Pow(Math.Sin(Math.PI / (2.0 * m)), 2);
            Max = -4 * h2 * Math.Pow(Math.Cos(Math.PI / (2.0 * n)), 2) - 4 * k2 * Math.Pow(Math.Cos(Math.PI / (2.0 * m)), 2);
            return 2 / (Max + Min);
        }
    }
}
