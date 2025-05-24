using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Reflection.Emit;
using System.Text;
using System.Windows.Forms;
using static System.Windows.Forms.VisualStyles.VisualStyleElement.Button;

namespace LAB_CHM_2023_3_1
{
    public partial class Form1 : Form
    {
        public Form1()
        {
            InitializeComponent();
        }

        private void тестоваяЗадачаToolStripMenuItem_Click(object sender, EventArgs e)
        {
            Form3 form3 = new Form3();
            form3.Show();
        }

        private void основнаяЗадачаToolStripMenuItem_Click(object sender, EventArgs e)
        {
            Form2 form2 = new Form2();
            form2.Show();
        }

        private void label46_Click(object sender, EventArgs e)
        {
        }

        private void groupBox2_Enter(object sender, EventArgs e)
        {
        }

        private void button1_Click(object sender, EventArgs e)
        {
            int n = Convert.ToInt32(textBox1.Text);
            int m = Convert.ToInt32(textBox2.Text);
            int N_max = Convert.ToInt32(textBox3.Text);
            double Eps = Convert.ToDouble(textBox4.Text);
            double Tau = Convert.ToDouble(textBox28.Text);

            Solution solution = new Solution();
            var result = solution.SolveTestTask(n, m, N_max, Eps, Tau);

            // Обновление массивов для передачи в отрисовку
            v1 = result.V1;
            u = result.U;

            // Заполнение таблиц
            FillDataGridView(dataGridView1, dataGridView2, dataGridView3, n, m, result.X, result.Y, u, v1, result.MaxPogr, result.XMax, result.YMax);

            // Справка
            textBox9.Text = Convert.ToString(result.P);
            textBox10.Text = Convert.ToString(result.EpsMax);
            textBox11.Text = Convert.ToString(result.MaxPogr);
            textBox15.Text = Convert.ToString(result.MaxF);
            textBox16.Text = Convert.ToString(result.MaxR1);
            textBox12.Text = Convert.ToString(result.XMax);
            textBox13.Text = Convert.ToString(result.YMax);
            textBox14.Text = "Нулевое начальное приближение";
        }

        private void button2_Click(object sender, EventArgs e)
        {
            int n = Convert.ToInt32(textBox8.Text);
            int m = Convert.ToInt32(textBox7.Text);
            int N_max = Convert.ToInt32(textBox6.Text);
            double Eps = Convert.ToDouble(textBox5.Text);
            double Tau = Convert.ToDouble(textBox29.Text);

            Solution solution = new Solution();
            var result = solution.SolveMainTask(n, m, N_max, Eps, Tau);

            // Обновление массивов для передачи в отрисовку
            v2 = result.V2;
            v2_2 = result.V2_2;

            // Заполнение таблиц
            FillDataGridView(dataGridView4, dataGridView5, dataGridView6, n, m, result.X, result.Y, v2, v2_2, result.MaxPogr, result.XMax, result.YMax);

            // Справка
            textBox17.Text = Convert.ToString(result.P);
            textBox18.Text = Convert.ToString(result.EpsMax);
            textBox19.Text = Convert.ToString(result.MaxF);
            textBox20.Text = Convert.ToString(result.MaxR1);
            textBox26.Text = Convert.ToString(result.P2);
            textBox27.Text = Convert.ToString(result.EpsMax2);
            textBox21.Text = Convert.ToString(result.MaxF2);
            textBox22.Text = Convert.ToString(result.MaxR);
            textBox23.Text = Convert.ToString(result.MaxPogr);
            textBox24.Text = Convert.ToString(result.XMax);
            textBox25.Text = Convert.ToString(result.YMax);
        }

        private void button3_Click_1(object sender, EventArgs e)
        {
            int n = Convert.ToInt32(textBox1.Text);
            int m = Convert.ToInt32(textBox2.Text);
            Form4 form4 = new Form4(u, v1, n, m, false);
            form4.Show();
        }

        private void button4_Click(object sender, EventArgs e)
        {
            int n = Convert.ToInt32(textBox8.Text);
            int m = Convert.ToInt32(textBox7.Text);
            Form4 form4 = new Form4(v2, v2_2, n, m, true);
            form4.Show();
        }

        private void button5_Click(object sender, EventArgs e)
        {
            int n = Convert.ToInt32(textBox1.Text);
            int m = Convert.ToInt32(textBox2.Text);

            Solution solution = new Solution();
            double tau = solution.CalculateTau(n, m);

            MessageBox.Show(
                "Оптимальное T(тау)  = " + tau.ToString(),
                "Подсчёт ТАУ");
            textBox28.Text = Convert.ToString(tau);
        }

        private void button6_Click(object sender, EventArgs e)
        {
            int n = Convert.ToInt32(textBox8.Text);
            int m = Convert.ToInt32(textBox7.Text);

            Solution solution = new Solution();
            double tau = solution.CalculateTau(n, m);

            MessageBox.Show(
                "Оптимальное T(тау)  = " + tau.ToString(),
                "Подсчёт ТАУ");
            textBox29.Text = Convert.ToString(tau);
        }

        private void FillDataGridView(DataGridView dataGridView1, DataGridView dataGridView2, DataGridView dataGridView3, int n, int m, double[] x, double[] y, double[][] u, double[][] v, double MaxPogr, double xMax, double yMax)
        {
            dataGridView1.Rows.Clear();
            dataGridView1.Columns.Clear();
            dataGridView1.Columns.Add("C1", "");
            dataGridView1.Columns[0].Width = 50;
            dataGridView1.Columns[0].Frozen = true;
            dataGridView1.Columns.Add("C2", "i");
            dataGridView1.Columns[1].Width = 50;
            dataGridView1.Columns[1].Frozen = true;

            dataGridView2.Rows.Clear();
            dataGridView2.Columns.Clear();
            dataGridView2.Columns.Add("C2", "");
            dataGridView2.Columns[0].Width = 50;
            dataGridView2.Columns[0].Frozen = true;
            dataGridView2.Columns.Add("C3", "i");
            dataGridView2.Columns[1].Width = 50;
            dataGridView2.Columns[1].Frozen = true;

            dataGridView3.Rows.Clear();
            dataGridView3.Columns.Clear();
            dataGridView3.Columns.Add("C4", "");
            dataGridView3.Columns[0].Width = 50;
            dataGridView3.Columns[0].Frozen = true;
            dataGridView3.Columns.Add("C5", "i");
            dataGridView3.Columns[1].Width = 50;
            dataGridView3.Columns[1].Frozen = true;

            for (int i = 0; i <= n; i++)
            {
                dataGridView1.Columns.Add(Convert.ToString(i), Convert.ToString(i));
                dataGridView2.Columns.Add(Convert.ToString(i), Convert.ToString(i));
                dataGridView3.Columns.Add(Convert.ToString(i), Convert.ToString(i));
            }

            dataGridView1.Rows.Add("j", "Y\\X");
            dataGridView2.Rows.Add("j", "Y\\X");
            dataGridView3.Rows.Add("j", "Y\\X");

            for (int i = 0; i <= n; i++)
            {
                dataGridView1.Columns[i + 2].HeaderText = i.ToString();
                dataGridView2.Columns[i + 2].HeaderText = i.ToString();
                dataGridView3.Columns[i + 2].HeaderText = i.ToString();

                dataGridView1.Rows[0].Cells[i + 2].Value = x[i];
                dataGridView2.Rows[0].Cells[i + 2].Value = x[i];
                dataGridView3.Rows[0].Cells[i + 2].Value = x[i];
            }

            for (int j = 0; j <= m; j++)
            {
                dataGridView1.Rows.Add();
                dataGridView2.Rows.Add();
                dataGridView3.Rows.Add();

                for (int i = 0; i <= 1; i++)
                {
                    dataGridView1.Rows[j + 1].Cells[0].Value = j;
                    dataGridView1.Rows[j + 1].Cells[1].Value = y[j];
                    dataGridView2.Rows[j + 1].Cells[0].Value = j;
                    dataGridView2.Rows[j + 1].Cells[1].Value = y[j];
                    dataGridView3.Rows[j + 1].Cells[0].Value = j;
                    dataGridView3.Rows[j + 1].Cells[1].Value = y[j];
                }
            }

            double Pogr;
            for (int j = 0; j <= m; j++)
            {
                for (int i = 0; i <= n; i++)
                {
                    Pogr = Math.Abs(u[i][j] - v[i][j]);
                    v[i][j] = Math.Round(v[i][j] * 1000) / 1000;
                    u[i][j] = Math.Round(u[i][j] * 1000) / 1000;

                    dataGridView1.Rows[j + 1].Cells[i + 2].Value = u[i][j];
                    dataGridView2.Rows[j + 1].Cells[i + 2].Value = v[i][j];
                    dataGridView3.Rows[j + 1].Cells[i + 2].Value = Pogr;
                }
            }
        }

        // Объявление массивов для передачи в отрисовку
        double[][] v1;
        double[][] u;
        double[][] v2;
        double[][] v2_2;
    }
}
