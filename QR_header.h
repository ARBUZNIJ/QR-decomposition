﻿#pragma once
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <omp.h>
#include <chrono>
//#include <C:\Program Files (x86)\Intel\oneAPI\mkl\2023.2.0\include\mkl.h>

//const size_t bs = 32;
using namespace std;
template <typename T>
class QR
{
private:

	size_t i, j, k, n, m, block_size;
	T* Q, * R, * REF_A, * v, * u, * factor_block, * factor, * w, * z, * R_tmp;
	T eps = 1e-10;

public:

	QR(bool flag, size_t size, size_t b_s) : n(size), block_size(b_s) //flag = 0 - ââîä ðàíäîìíûõ ÷èñåë, èíà÷å ââîä ñ êëàâèàòóðû
	{
		Q = new T[n * n]();
		R = new T[n * n];
		R_tmp = new T[n * n];
		REF_A = new T[n * n];
		u = new T[n];
		v = new T[n * block_size];
		factor_block = new T[n];
		factor = new T[n * n];
		w = new T[n * n](); //ПРОМАЗАЛИ С РАЗМЕРОМ
		z = new T[n];

		if (flag)
		{

			cout << "Enter matrix:" << endl;
			for (i = 0; i < n; i++)
				for (j = 0; j < n; j++)
					cin >> R[i * n + j];

			for (i = 0; i < n; i++)

				copy(R + i * n, R + (i + 1) * n, REF_A + i * n);

		}
		else {


			for (i = 0; i < n; i++)
				for (j = 0; j < n; j++)

					R[i * n + j] = T(rand()) / RAND_MAX;

			for (i = 0; i < n; i++)
				R[i * n + i] += n;

			//#pragma omp parallel for schedule(static) private(i)
			for (i = 0; i < n; i++)
				copy(R + i * n, R + (i + 1) * n, REF_A + i * n);

			/*auto start{ chrono::steady_clock::now() };
			LAPACKE_dgeqrf(LAPACK_ROW_MAJOR, n, n, REF_A, n, v);
			LAPACKE_dorgqr(LAPACK_ROW_MAJOR, n, n, n, REF_A, n, v);
			auto end{ chrono::steady_clock::now() };
			chrono::duration<double> elapsed_seconds = end - start;
			cout << "Time spent (MKL): " << elapsed_seconds.count() << endl;*/

			for (i = 0; i < n; i++)
				copy(R + i * n, R + (i + 1) * n, REF_A + i * n);

		}

	}

	double* returnR()
	{
		return R;
	}

	double* returnV()
	{
		return v;
	}

	void count_v_gamma(size_t column, size_t num_in_block)
	{
		T scl, gamma;
		//#pragma omp parallel for simd reduction(+: scl) 
		scl = 0;

		for (size_t i = 0; i < n - column; i++)
		{
			u[i] = R[(i + column) * n + column];
			scl += u[i] * u[i];
		}

		if (scl < eps)
		{
			v[num_in_block * n + column] = 1;
			gamma = 0.5;
			return;
		}

		else
		{
			scl = 1 / sqrt(scl);
			u[0] *= scl;
			gamma = (1 + abs(u[0]));
			v[num_in_block * n + column] = sgn(u[0]) * gamma;

			//#pragma omp parallel for simd
			for (size_t i = column + 1; i < n; i++)
				v[num_in_block * n + i] = u[i - column] * scl;
		}
	}

	T sgn(T val)                //íåêàíîíè÷íûé sign(x), êîòîðûé âîçâðàùàåò 1 äëÿ íåîòðèöàòåëüíûõ x, èíà÷å -1
	{
		if (val >= 0)
			return 1;
		else return -1;
	}

	T scal(size_t v_num_in_block, size_t v_ind, size_t a_col)  //ôóíêöèÿ ñêàëÿðíîãî ïðîèçâåäåíèÿ äëÿ k ñòîëáöa ìàòðèöû A è âåêòîðà v, ïåðâûé íåíóëåâîé ýëåìåíò â êîòîðîì ðàñïîëàãàåòñÿ â ïîçèöèè ind
	{
		T res = 0;
		//#pragma omp simd
		for (size_t i = v_ind; i < n; i++)
			res += R[i * n + a_col] * v[v_num_in_block * n + i];
		return res;
	}
	void HHolder_R()
	{
		size_t m = 0;

		for (; m < (n / block_size); m++)

			HHolder_Block(m * block_size, block_size);

		if (n % block_size != 0)
			HHolder_Block(m * block_size, n - m * block_size);
	}
	void HHolder_Block(size_t i_start, size_t b_size)
	{
		size_t num_in_block;

		for (i = 0; i < block_size; i++)
			memset(v + i * n, 0, n * sizeof(T));

		for (j = i_start; j < i_start + b_size; j++)
		{
			num_in_block = j - i_start;
			count_v_gamma(j, num_in_block); //âû÷èñëèëè vi è çàïîëíèëè èìè v[]

			for (k = j; k < i_start + b_size; k++)
				factor_block[k - j] = scal(num_in_block, j, k) / abs(v[num_in_block * n + j]); //âû÷èñëèëè ê-òû äëÿ vi â áëîêå

			for (i = j; i < n; i++)
			{
				for (k = j; k < i_start + b_size; k++)
					R[i * n + k] -= v[num_in_block * n + i] * factor_block[k - j]; //ïîâû÷èòàëè èç áëîêà
			}
		}
		

		//for (i = 0; i < block_size; i++) ///?
		memset(w, 0, n * n * sizeof(T));

		double norm = 0.0;
		for (i = i_start; i < n; i++)
		{
			norm += v[i] * v[i];
		}
		norm = -2 / (norm);

		for (i = i_start; i < n; i++)
		{
			w[i * b_size] = (norm)*v[i];
		}

		for (i = 1; i < b_size; i++) //ïî ñòîëáöàì ðåçóëüòàòà â w
		{

			norm = 0.0;
			for (k = 0; k < n; k++)
			{
				norm += v[i * n + k] * v[i * n + k];
			}
			norm = -2 / (norm);

			for (j = 0; j < n; j++) //ïî ñòðîêàì W
			{
				for (k = i; k < n; k++) //ïî ñòðîêàì Y == v
				{
					for (m = 0; m < i; m++) //ïî ñòîëáöàì W
					{
						w[j * b_size + i] += (w[j * b_size + m] * v[m * n + k] + (j == k ? 1 : 0)) * v[i * n + k] * (norm);
					}
				}
			}
		}

		//for (i = 0; i < n; i++)
		memset(Q, 0, n * n * sizeof(T));

		for (j = 0; j < n; j++) //ïî ñòðîêàì W
		{
			for (k = 0; k < n; k++) //ïî ñòðîêàì Y == v
			{
				for (m = 0; m < b_size; m++) //ïî ñòîëáöàì W
				{
					Q[j * n + k] += v[m * n + j] * w[k * b_size + m];
				}
			}
		}

		for (j = 0; j < n; j++)
			Q[j * n + j] += 1;

		copy(R, R + n * n, R_tmp);
		for (i = 0; i < n; i++)
			memset(R + i * n + i_start + b_size, 0, (n - i_start - b_size) * sizeof(T)); //до конца строки всё занулили 

		for (j = 0; j < n; j++)
		{
			for (m = 0; m < n; m++)
			{
				for (k = i_start + b_size; k < n; k++)
				{
					R[j * n + k] += Q[j * n + m] * R_tmp[m * n + k];
				}
			}
		}
	}

	void HHolder_Q()
	{
		//#pragma omp parallel for private(i)
		for (i = 0; i < n; i++)
			copy(REF_A + i * n, REF_A + (i + 1) * n, Q + i * n);

		for (i = 0; i < n; i++)
		{
			//#pragma omp parallel for private(j)
			for (j = 0; j < n; j++)
				Q[j * n + i] /= R[i * n + i];

			copy(R + i * (n + 1) + 1, R + (i + 1) * n, factor);

			//#pragma omp parallel for private(k)
			for (k = 0; k < n; k++)
				for (j = i + 1; j < n; j++)

					Q[k * n + j] -= factor[j - (i + 1)] * Q[k * n + i];

		}
	}


	bool check()
	{
		T* tmp = new T[n * n]();

		for (i = 0; i < n; i++)
			for (k = 0; k < n; k++)
				for (j = 0; j < n; j++)
					tmp[i * n + j] += Q[i * n + k] * R[k * n + j];

		for (i = 0; i < n; i++)
			for (j = 0; j < n; j++)
				if (abs(REF_A[i * n + j] - tmp[i * n + j]) >= eps)
					return false;

		return true;
	}

	void out(char s)
	{
		cout << endl;
		if (s == 'R')

			for (i = 0; i < n; i++)
			{
				for (j = 0; j < n; j++)
					cout << R[i * n + j] << ' ';
				cout << endl;

			}

		else  if (s == 'Q')

			for (i = 0; i < n; i++)
			{
				for (j = 0; j < n; j++)
					cout << Q[i * n + j] << ' ';
				cout << endl;
			}
		cout << endl;
	}

	void transpQ()
	{
		for (i = 0; i < n; i++)

			//#pragma omp parallel for private(j)
			for (j = i + 1; j < n; j++)
				swap(Q[i * n + j], Q[j * n + i]);

		return;
	}

	~QR()
	{
		delete[]Q;
		delete[]R;
		delete[]v;
		delete[]u;
		delete[]REF_A;
		delete[]factor;
		delete[]factor_block;
		delete[]w;
		delete[]z;
		delete[]R_tmp;

		Q = R = u = REF_A = v = factor = nullptr;
		i = j = k = m = n = block_size = 0;
	}
};