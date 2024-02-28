#include <omp.h>
#include "QR_header.h"

int main()
{
	//int count;
	size_t n, block_size;
	cout << "Enter matrix dimentions: " << endl;
	cin >> n;
	/*while (1)
	{*/
		cout << "Enter block size: ";
		cin >> block_size;
		//count = 0;
		/*while (count < 5)
		{*/

			QR <double> t(1, n, block_size);
			/*auto start_A{ chrono::steady_clock::now() };*/

			t.HHolder_R();

			t.out('R');

			/*auto end_A{ chrono::steady_clock::now() };
			chrono::duration<double> elapsed_seconds_A = end_A - start_A;
			auto start_Q = chrono::steady_clock::now();*/

			t.HHolder_Q();

			t.out('Q');

			/*auto end_Q = chrono::steady_clock::now();
			chrono::duration<double>  elapsed_seconds_Q = end_Q - start_Q;
			cout << "Time spent: " << elapsed_seconds_A.count() << " sec + " << elapsed_seconds_Q.count() << " sec = ";
			elapsed_seconds_A += elapsed_seconds_Q;
			cout << elapsed_seconds_A.count() << endl;
			count++;*/

			if (t.check())

				cout << "QR-decomposition is correct";
			else cout << "QR-decomposition is incorrect";

			cout << endl;
		/*}*/
	//}
	/*t.HHolder_Q();
	end = chrono::steady_clock::now();
	elapsed_seconds = end - start;
	cout << "Time spent: " << elapsed_seconds.count() << " sec" << endl;
	if (t.check())
		cout << "QR-decomposition is correct";
	else cout << "QR-decomposition is incorrect";
	cout << endl;*/
	/*QR <double> t(0, n, block_size);
	t.out('R');*/

			/*void norm_v()
			{
				memset(factor_block, 0, block_size * sizeof(T));

				for (i = 0; i < block_size; i++)
					for (j = 0; j < n; j++)
						factor_block[i] += v[i * n + j] * v[i * n + j];

				for (i = 0; i < block_size; i++)
					factor_block[i] = (-2) / factor_block[i];

				for (i = 0; i < block_size; i++)
					for (j = 0; j < n; j++)
						v[i * n + j] *= factor_block[i];
			}*/

	return 0;
}