using System;
using System.Text;
using static System.Runtime.InteropServices.JavaScript.JSType;
namespace DataAnalysisLab1;

static class Program
{
    /// <summary>
    ///  Nhe main entry point for the application.
    /// </summary>
    [STAThread]
    static void Main()
    {
        Console.InputEncoding = Encoding.Unicode;
        Console.OutputEncoding = Encoding.Unicode;
        double q = 6;
        double M = Math.Pow(10, q);
        double b = 37;
        int N = 170;
        double a0 = 3;

        double[][] period_Array = new double[(int)(N)][];

        for (int i = 0; i < N; i++)
        {
            period_Array[i] = Left_Metod(a0, b, M, N);
            a0 = period_Array[i][N - 1] * M;
        }

        //X_Square(N, 21, period_Array[0]);
        double[] NormalArray = new double[(int)(N)];
        for (int i = 0; i < N; i++)
        {
            Console.Write($"{i + 1}-е нормально роподілене число: ");
            NormalArray[i] = ModelingNormalDistribution(N, period_Array[i]);
        }

        //Dictionary<string, double> NumberParams = new Dictionary<string, double>();

        double x_max, x_min;
        double[] minmax = MaxMin(NormalArray);
        x_min = minmax[0];
        x_max = minmax[1];
        double R_range = x_max - x_min;
        int k_intervals = 17;
        double h_width = R_range / k_intervals;

        double[] x_intervalLimits = IntervalLimits_Find(k_intervals, h_width, x_max, x_min);
        double[] n_Rates = IntervalRates(k_intervals, N, NormalArray, x_min, h_width);
        double[] x_discVar = DiscreteVarianta(k_intervals, h_width, x_intervalLimits);

        double x_selAve = SelectiveAverage(N, k_intervals, n_Rates, x_discVar);
        double sigma = SigmaCalc(N, k_intervals, x_selAve, x_discVar, n_Rates);
        double selDisp_fixed = Fixed_SelectiveDispersion(N, k_intervals, x_selAve, x_intervalLimits);

        ThreeSigma_Rule(sigma, x_selAve, NormalArray, N);

        double[] empiric = EmpiricDistribution(n_Rates, k_intervals, N);

        //double[] mi_TeoreticalRates_Find = new double[k_intervals];
        //for (int i = 0; i < k_intervals-1; i++)
        //{
        //    mi_TeoreticalRates_Find[i] = TeoreticalRates_Find(x_intervalLimits[i], x_intervalLimits[i + 1], N);
        //}
        for (int i = 0; i < k_intervals-1; i++)
        {
            double step = 4.0 / k_intervals;
            Console.Write($"Interval start: {(2.0 + step * (i)):F16} | ");
            TeoreticalRates_Find(2.0 + step * (i),
                2.0 + step * (i + 1), N);
        }


        //YastremskyiCriteriy(empiric, mi_TeoreticalRates_Find, k_intervals, N);
    }

    public static double[] Left_Metod(double a0, double b, double M, int N)
    {
        double[] period_Array = new double[(int)(N)];
        period_Array[0] = (a0 / M);//ready of random numbers

        double[] ai_Array = new double[(int)(N)];
        ai_Array[0] = a0;

        for (int i = 1; i < N; i++)
        {
            ai_Array[i] = (b * ai_Array[i - 1]) % M;
            period_Array[i] = (ai_Array[i] / M);
            //Console.WriteLine(period_Array[i]);
        }

        //var dict = new Dictionary<double, int>();
        //foreach (var value in period_Array)
        //{
        //    if (dict.ContainsKey(value))
        //    {
        //        dict[value]++;
        //    }
        //    else
        //    {
        //        dict[value] = 1;
        //    } 
        //}
        //Console.WriteLine(dict.Count);

        return period_Array;
    }

    public static void X_Square(int N, int k, double[] test_Array)
    {
        int[] intervalSizes = new int[k];
        int interval;
        double p = 1.0 / k;//size of range for each interval
        //Elements in each interval
        for (int i = 0; i < N; i++)
        {
            //by deviding by size of intervals we get the interval which it belongs to
            interval = ((int)(test_Array[i] / p));
            intervalSizes[interval] += 1;
        }

        double X_Exp = 0;
        for (int i = 0; i < k; i++)
        {
            X_Exp += (
                ((intervalSizes[i] - p * N) * (intervalSizes[i] - p * N)) /
                    (p * N));
        }

        double X_teor = 31.4;
        if (X_Exp > X_teor)
        {
            Console.WriteLine($"{X_Exp}>{X_teor}");
            Console.WriteLine("Генератор неякісний");
        }
        else
        {
            Console.WriteLine($"{X_Exp}<{X_teor}");
            Console.WriteLine("Генератор якісний");
        }
    }

    public static double ModelingNormalDistribution(double N, double[] Numbers)
    {
        double ksi = 0;
        for (int i = 0; i < N; i++)
        {
            ksi += (Numbers[i]);
        }
        ksi -= N / 2;
        ksi = ksi * Math.Sqrt(12.0 / N);


        double D = 0.7;

        double a = 7;
        double nu = 0;
        nu = a + Math.Sqrt(D) * ksi;
        Console.WriteLine($"Число мю: {nu}");
        return nu;
    }


    public static double[] MaxMin(double[] distribution)
    {
        double x_min = distribution[0];
        double x_max = distribution[0];
        for (int i = 1; i < distribution.Length; i++)
        {
            if (distribution[i] < x_min)
            {
                x_min = distribution[i];
            }
            else if (distribution[i] > x_max)
            {
                x_max = distribution[i];
            }
        }
        double[] result = [x_min, x_max];

        return result;
    }

    public static double[] IntervalLimits_Find(int k_intervals, double h_width, double x_max, double x_min)
    {

        double[] x_intervalLimits = new double[k_intervals];
        for (int i = 0; i < x_intervalLimits.Length; i++)
        {
            x_intervalLimits[i] = x_min + (i) * h_width;
        }

        return x_intervalLimits;
    }

    public static double[] IntervalRates(int k_intervals, double N, double[] distribution, double x_min, double h_width)
    {
        double[] n_Rates = new double[k_intervals];
        for (int i = 0, interval; i < N; i++)
        {
            //by deviding by size of intervals we get the interval which it belongs to
            interval = ((int)((distribution[i] - x_min) / h_width - 0.0001));//0.0001 is needed because when time comes to max element we get the value equal to our number of intervals, which cause "index was outside the bounds of the array" error 
            n_Rates[interval] += 1;
        }

        Console.WriteLine("Частоти: ");
        for (int interval = 0; interval < k_intervals; interval++)
        {
            Console.WriteLine($"{interval + 1}: {n_Rates[interval]}");
        }

        return n_Rates;
    }


    public static double[] DiscreteVarianta(double k_intervals, double h_width, double[] x_intervalLimits)
    {
        double[] x_discVar = new double[(int)k_intervals];//discrete series varianta
        for (int i = 0; i < k_intervals; i++)
        {
            x_discVar[i] = x_intervalLimits[i] + h_width / 2;
        }

        return x_discVar;
    }

    public static double SelectiveAverage(double N, double k_intervals, double[] n_Rates, double[] x_discVar)
    {
        double x_selAve = 0;
        for (int i = 0; i < k_intervals; i++)
        {
            x_selAve += x_discVar[i] * n_Rates[i];
        }
        x_selAve /= N;//selective average
        Console.WriteLine($"Середнє вибіркове: {x_selAve}");

        return x_selAve;
    }

    public static double SigmaCalc(double N, double k_intervals, double x_selAve, double[] x_discVar, double[] n_Rates)
    {
        double sigmaSqv = 0;
        for (int i = 0; i < k_intervals; i++)
        {
            sigmaSqv += (x_discVar[i] - x_selAve) * (x_discVar[i] - x_selAve) * n_Rates[i];
        }
        sigmaSqv /= N;

        double sigma = Math.Sqrt(sigmaSqv);//average square deviation
        Console.WriteLine($"Середньоквадратичне відхилення: {sigma}");

        return sigma;
    }

    public static double Fixed_SelectiveDispersion(double N, int k_intervals, double x_selAve, double[] x_intervalLimits)
    {
        double selDisp_fixed = 1 / (N - 1);
        for (int i = 0; i < k_intervals; i++)
        {
            selDisp_fixed += (x_intervalLimits[i] - x_selAve) * (x_intervalLimits[i] - x_selAve);
        }

        selDisp_fixed /= (N - 1);
        Console.WriteLine($"Виправлена вибіркова дисперсія: {selDisp_fixed}");

        return selDisp_fixed;
    }


    public static void ThreeSigma_Rule(double sigma, double x_Ave, double[] distribution, int N)
    {
        double alpha = 0.05;
        for (int i = 0; i < N; i++)
        {
            if ((Math.Abs(distribution[i] - x_Ave)) > 3 * sigma)
            {
                //Console.WriteLine("Правило трьох сигм порушено!");
                Console.WriteLine($"{i + 1}-й елемент розподілу має аномальне значення:{distribution[i]}");
            }
        }

        //Console.WriteLine("Правило трьох сигм не порушено!");
    }

    public static double[] EmpiricDistribution(double[] n_Rates, int k_intervals, double N)
    {
        double[] sizes = new double[k_intervals];
        double sum = 0;
        Console.WriteLine("Емпіричний розподіл: ");

        for (int i = 0; i < k_intervals; i++)
        {
            sum += n_Rates[i];
            sizes[i] = sum / N;
        }

        double heightInterval = 0.05;

        int interval = k_intervals - 1;

        Console.WriteLine(new string(' ', 8) + "^");
        for (double currentHeight = 1.0; currentHeight > 0; currentHeight -= heightInterval)
        {
            Console.Write(new string(' ', 3));
            if (sizes[interval] > currentHeight - heightInterval)
            {
                Console.Write($"{sizes[interval],-5:F2}" + "|");
                Console.WriteLine(new string(' ', (interval) * 3) + new string('_', 4));
                interval--;
            }
            else
            {
                Console.WriteLine(new string(' ', 5) + "|");
            }
        }
        Console.WriteLine(new string('_', 8) + "|" + new string('_', k_intervals * 3) + '>');
        Console.WriteLine(new string(' ', 8) + "|" + new string(' ', k_intervals * 3) + 'x');
        Console.WriteLine(new string(' ', 8) + "|");
        Console.WriteLine(new string(' ', 8) + "|");

        return sizes;
    }


    public static double TeoreticalRates_Find(double a_intevalStart, double b_intervalEnd, int N_total)
    {
        double Range = b_intervalEnd - a_intevalStart;
        double h_width = Range / N_total;
        double x0 = a_intevalStart;

        double tempSum = 0;
        double x = x0 + h_width / 2;
        double xi;

        //double derivative;
        //double max_xiDerive = double.MinValue;
        for (double i = 0; i < N_total - 1; i++)
        {
            xi = x + h_width * i;
            tempSum += (LaplasFunction(xi));

            //derivative = Math.Abs(SecondLaplasDerivative(xi));
            //if (max_xiDerive < derivative)
            //{
            //    max_xiDerive = derivative;
            //}
        }

        //double R_error = ((b_intervalEnd - a_intevalStart) / 24) * max_xiDerive * Math.Pow(h_width, 2);
        double teoreticalRate = tempSum * h_width * N_total;

        //Console.Write($"Найбільша похідна: {max_xiDerive} | ");
        Console.WriteLine($"Теоретична частота: {teoreticalRate}");
        return teoreticalRate;
    }

    public static void YastremskyiCriteriy(double[] n_empiricRates, double[] m_teoreticRates, int k_intervals, int N_total)
    {
        double Q = 0;
        for (int i = 1; i < k_intervals; i++)
        {
            Q += (Math.Pow(n_empiricRates[i] - m_teoreticRates[i], 2) /
                    (m_teoreticRates[i] * (1 - m_teoreticRates[i] / N_total))
                    );
        }

        double beta = 0.6;
        double J = Math.Abs(Q - k_intervals) /
            Math.Sqrt(2 * k_intervals + 4 * beta);

        if (J <= 3) Console.WriteLine("\nГіпотеза приймається");
        else Console.WriteLine("\nГіпотеза не приймається");
    }

    public static double LaplasFunction(double x)
    {
        //double delta = 0.0001;
        double result = (Math.Exp((
                            -(x * x) / 2))
                        );
        return result;
    }

    public static double SecondLaplasDerivative(double x)
    {
        return (
                (Math.Pow(x, 2) - 1) * Math.Exp(-(x * x) / 2)
                                                            );
    }
}

