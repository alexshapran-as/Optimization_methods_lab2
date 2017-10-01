#include <iostream>
#include <fstream>
#include "simplex.h"
#include <vector>
#include <cstdlib>
#include <algorithm>
#include <string>
#include <iomanip>

void set_data(std::vector<double> &c, std::vector<double> &b, double** a, unsigned int &n, unsigned int &m)
{
	using std::cout;
	using std::endl;
	using std::vector;

	double temp = {0.0};
	std::ifstream fin("lab2.txt");
	if (!fin.is_open())
	{
		cout << endl << "Ошибка открытия файла.";
		exit(EXIT_FAILURE);
	}
	fin >> min_max;
	if (min_max != "min" && min_max != "max")
	{
		cout << endl << "Неправильно указано, к чему стремится ЦФ F.";
		exit(EXIT_FAILURE);
	}
	while (fin >> temp)
	{
		if (min_max == "min")
			c.push_back(temp);
		else if (min_max == "max")
			c.push_back(-temp);
		n++;
	}
	fin.clear();
	while (fin.get() != '\n');
	while (fin >> temp)
	{
		b.push_back(temp);
		m++;
	}
	fin.clear();
	while (fin.get() != '\n');

	if (m > 2*n)
	{
		cout << endl << "Число линейно независимых уравнений больше числа переменных."
		" Такая СЛАУ несовместна.";
		exit(EXIT_FAILURE);
	}

	for (unsigned int i = {0}; i < m; ++i)
		a[i] = new double [n];
	for (unsigned int i = {0}; i < m; ++i)
	{
		for (unsigned int j = {0}; j < n; ++j)
		{
			if (!(fin >> a[i][j]))
			{
				cout << endl << "Не удалось считать данные для вектора A = [a1,..,an]!";
				exit(EXIT_FAILURE);
			}
		}
	}
	fin.close();
	show_task(c,b,a);
}

void show_task(std::vector<double> c, std::vector<double> b, double** a)
{
	using std::cout;
	using std::cin;
	using std::endl;

	cout << endl << "Постановка задачи ЛП. Требуется найти решение следующей задачи:" 
	<< endl << std::setw(26) << "F = cx -> " <<  min_max << endl 
	<< std::setw(25) << "Ax <= b" << endl << std::setw(25) << "x >= 0" << endl
	<< "Каноническая форма задачи ЛП:" << endl;
	if (min_max == "min")
		cout << std::setw(17) << "F = ";
	else
		cout << std::setw(17) << "-F = ";
	for (unsigned int i = {0}; i < c.size(); ++i)
	{
		if (c[i] < 0)
			cout << c[i] << "x" << i+1;
		else
			cout << "+" << c[i] << "x" << i+1;
	}
	cout << " -> ";
	if (min_max == "min")
		cout <<  min_max << endl << endl; 
	else
		cout <<  "min" << endl << endl; 
	unsigned int k = c.size() + 1;
	for (unsigned int i = {0}; i < b.size(); ++i)
	{
		cout << "                 ";
		for (unsigned int j = {0}; j < c.size(); ++j)
		{
			if (a[i][j] != 0)
			{
				if (a[i][j] != 1 && a[i][j] != -1)
				{
					if (a[i][j] < 0 || j == 0)
						cout << a[i][j] << "x" << j+1;
					else
						if (a[i][j-1] != 0)
							cout << "+" << a[i][j] << "x" << j+1;
						else 
							cout << a[i][j] << "x" << j+1;
				}
				else
				{
					if (a[i][j] < 0)
						cout << "-x" << j+1;
					else if (j != 0)
						cout << "+" << "x" << j+1;
					else 
						cout << "x" << j+1;
				}
		
			}		
		}
		cout << "+x" << k++ << " = " << b[i] << endl;
	}
	cout << endl << "             ";
	for (unsigned int i = {0}; i < c.size()+b.size(); ++i)
	{
		if (i != c.size()+b.size()-1)
			cout << "x" << i+1 << ", ";
		else
			cout << "x" << i+1;
	}
	cout << " >= 0" << endl << "Алгоритм преобразования симплекс-таблицы (жордановы исключения):";

}

double * simplex_method(std::vector<double> c, std::vector<double> b, double** a)
{
	using std::cout;
	using std::endl;
	using std::vector;
	using std::string;

	vector<string> x(b.size() + c.size() + 1);

	double ** arr = new double * [b.size()+1];  // Двумерный массив - симплекс матрица
	for (unsigned int i = {0}; i < b.size()+1; ++i)
		arr[i] = new double [c.size()+1];
	for (unsigned int i = {0}; i < b.size(); ++i)
		arr[i][0] = b[i]; // Заполняем первый столбец свободных членов
	arr[b.size()][0] = {0}; // Свободный элемент в функции F
	for (unsigned int i = {0}; i < b.size(); ++i)
	{
		unsigned int k = {1};
		for (unsigned int j = {0}; j < c.size(); ++j)
		{
			arr[i][k] = a[i][j]; // Заполняем все остальные элементы, начиная со второго столбца, без последней строки
			k++;
		}
	} 
	unsigned int k = {1};
	for (unsigned  int i = {0}; i < c.size(); ++i)
	{
		arr[b.size()][k] = -c[i]; // Заполняем последнюю строку со второго столбца коэфф-ми при x в ф-ии F
		k++;
	}

	x = print(arr, b.size(), c.size());

	while (not_reference(arr,b.size()))
	{
		cout << endl << "Недопустимое решение. Производим замену базиса," 
		" чтобы получить опорное решение.";
		unsigned int i = {0};
		while (arr[i][0] > 0 && i < b.size() + 1)
			i++;
		unsigned int j = {1}; // Разрешающий столбец
		unsigned int k = {0};
		for (; j < c.size() + 1; ++j)
		{
			if (arr[i][j] < 0)
			{
				k++;
				break;
			}
		}
		if (k == 0)
		{
			cout << endl << "Задача не имеет допустимых решений!";
			exit(EXIT_FAILURE);
		}
		vector<double> temp;
		for (unsigned int i = {0}; i < b.size(); ++i)
		{
			if (arr[i][0] != 0 || arr[i][j] > 0)
				temp.push_back(arr[i][0] / arr[i][j]); // Заполняем временный вектор temp отношениями свободных членов к элементам разрешающего столбца
		}
		std::sort(temp.begin(),temp.end()); // Сортируем эти отношения
		unsigned int m = {0};
		while (temp[m] < 0)
			m++; // Ищем номер первого положительного элемента в отсортированном векторе temp 
		i = {0}; // Разрешающая строка
		if (arr[ri][rj] == 0)
		{
			cout << endl << "Деление на 0 недопустимо!";
			exit(EXIT_FAILURE);
		}
		for (; temp[m] != arr[i][0] / arr[i][j]; ++i);
		for (unsigned int m = {0}; m < b.size() + 1; ++m)
		{
			for (unsigned int n = {0}; n < c.size() + 1; ++n)
			{
				if (m != i && n != j) // Преобразование эл-та, который не стоит ни в разрешающей строке, ни в разрешающем столбце
				{
					arr[m][n] = arr[m][n] - (arr[m][j] * arr[i][n]) / (arr[i][j]);
				}
			}
		}
		for (unsigned int m = {0}; m < b.size() + 1; ++m)
		{
			for (unsigned int n = {0}; n < c.size() + 1; ++n)
			{
				if (m != i && n == j) // Преобразование эл-та, который стоит в разрешающем столбце, но не в разрешающей строке
				{
					if (arr[m][n] != 0)
						arr[m][n] = -arr[m][n] / arr[i][j];
				}
				if (m == i && n != j) // Преобразование эл-тва, который стоит в разрешающей строке, но не в разрешающем столбце
				{
					arr[m][n] = arr[m][n] / arr[i][j];
				}
			}
		}
		arr[i][j] = 1.0 / arr[i][j]; // Преобразование разрешающего эл-та
		ri = i; 
		rj = j;
		x = print(arr, b.size(), c.size());
	}

	cout << endl << "Найдено опорное решение. Проводим проверку на оптимальность.";
	while (not_optimal(arr,b.size(),c.size()))
	{
		cout << endl << "Неоптимальное решение. Производим замену базиса.";
		unsigned int j = {1}; // Новый разрешающий столбец 
		for (; j < c.size() + 1; ++j)
		{
			if (arr[b.size()][j] > 0) // Ищем разрешающий столбец
				break;
		}
		vector<double> temp;
		for (unsigned int i = {0}; i < b.size(); ++i)
		{
			if (arr[i][0] != 0 || arr[i][j] > 0)
				temp.push_back(arr[i][0] / arr[i][j]); 
		}
		std::sort(temp.begin(),temp.end());
		unsigned int m = {0};
		while (temp[m] < 0) // Аналогично ищем номер первого положительного отношения
			m++;
		unsigned int i = {0}; // Новая разрешающая строка
		if (arr[ri][rj] == 0)
		{
			cout << endl << "Деление на 0 недопустимо!";
			exit(EXIT_FAILURE);
		}
		for (; temp[m] != arr[i][0] / arr[i][j]; ++i);
		for (unsigned int m = {0}; m < b.size() + 1; ++m)
		{
			for (unsigned int n = {0}; n < c.size() + 1; ++n)
			{
				if (m != i && n != j)
				{
					arr[m][n] = arr[m][n] - (arr[m][j] * arr[i][n]) / (arr[i][j]);
				}
			}
		}
		for (unsigned int m = {0}; m < b.size() + 1; ++m)
		{
			for (unsigned int n = {0}; n < c.size() + 1; ++n)
			{
				if (m != i && n == j)
				{
					if (arr[m][n] != 0)
						arr[m][n] = -arr[m][n] / arr[i][j];
				}
				if (m == i && n != j)
				{
					arr[m][n] = arr[m][n] / arr[i][j];
				}
			}
		}
		arr[i][j] = 1.0 / arr[i][j];
		ri = i;
		rj = j;
		x = print(arr, b.size(), c.size()); // Вектор строк, в котором первыми идут свободные переменные, а затем идут базисные переменные
	}
	cout << endl << "Найдено оптимальное решение!" << endl;
	for (unsigned int i = {0}; i < c.size(); ++i)
		cout << x[i] << " = 0; ";
	unsigned int l = {0};
	for (unsigned int j = c.size(); j < b.size()+c.size(); ++j)
	{
		cout << endl << x[j] << " = " << arr[l++][0] << ";";
	}
	cout << endl;
	double result = {0.0};
	double sum = {0.0};
	if (min_max == "min") 
	{
		if (x[0].find("X") != -1)
		{
			cout << "F = " << arr[b.size()][0];
			result = arr[b.size()][0];
			cout << endl << endl << std::setw(35) << "Проверка:" << endl << endl << "       ";
			for (unsigned int i = {0}; i < b.size(); ++i)
			{
				for (unsigned int j = {0}; j < c.size(); ++j)
				{
					cout << arr[j][0] << "*" << a[i][j] << " + ";
					sum += arr[j][0] * a[i][j];
				}
				cout << " 0 = " <<  sum << endl << "       ";
				sum = {0.0};
			}
		}
		else
		{
			std::swap(arr[1][0],arr[2][0]);
			cout << "-Ф = " << -arr[b.size()][0];
			result = -arr[b.size()][0];
			cout << endl << endl << std::setw(35) << "Проверка:" << endl << endl << "       ";
			for (unsigned int i = {0}; i < b.size(); ++i)
			{
				for (unsigned int j = {0}; j < c.size(); ++j)
				{
					cout << arr[j][0] << "*" << a[i][j] << " + ";
					sum += arr[j][0] * a[i][j];
				}
				cout << " 0 = " <<  sum << endl << "       ";
				sum = {0.0};
			}
		}
	}
	else if (min_max == "max")
	{
		if (x[0].find("X") != -1)
		{
			cout << "F = " << -arr[b.size()][0]; // Инвертируем знак F, так как ищем максимум
			result = -arr[b.size()][0];
			cout << endl << endl << std::setw(35) << "Проверка:" << endl << endl << "       ";
			for (unsigned int i = {0}; i < b.size(); ++i)
			{
				for (unsigned int j = {0}; j < c.size(); ++j)
				{
					cout << arr[j][0] << "*" << a[i][j] << " + ";
					sum += arr[j][0] * a[i][j];
				}
				cout << " 0 = " <<  sum << endl << "       ";
				sum = {0.0};
			}
		}
		else
		{
			std::swap(arr[1][0],arr[2][0]);
			cout << "Ф = " << arr[b.size()][0];
			result = arr[b.size()][0];
			cout << endl << endl << std::setw(35) << "Проверка:" << endl << endl << "       ";
			for (unsigned int i = {0}; i < b.size(); ++i)
			{
				for (unsigned int j = {0}; j < c.size(); ++j)
				{
					cout << arr[j][0] << "*" << a[i][j] << " + ";
					sum += arr[j][0] * a[i][j];
				}
				cout << " 0 = " <<  sum << endl << "       ";
				sum = {0.0};
			}
		}
	}
	if (x[0].find("X") == -1)
	{
		cout << endl << endl << "Таким образом, решения" << endl << "прямой и двойственной задач совпадают:" 
		<< endl << "F = Ф = " << result;
	}


	for (unsigned int i = {0}; i < b.size() + 1; ++i)
		delete [] arr[i];
	delete [] arr;
}

bool not_reference(double ** arr, unsigned int m)
{
	unsigned int accept = {0};
	for (unsigned int i = {0}; i < m; ++i)
	{
		if (arr[i][0] < 0) // Проверка на допустимое решение
			accept++;
	}
	if (accept > 0)
		return true;
	else
		return false;
}

bool not_optimal(double ** arr, unsigned int m, unsigned int n)
{
	unsigned int check = {0};
	for (unsigned int j = {1}; j < n + 1; ++j)
	{
		if (arr[m][j] > 0) // Проверка на оптимальное решение
			check++;
	}
	if (check > 0)
		return true;
	else 
		return false;
}

std::vector<std::string> print(double ** arr, unsigned int m, unsigned int n)
{
	using std::cout;
	using std::endl;
	using std::vector;
	using std::string;

	cout << endl << "Симплекс-таблица:" << endl << std::setw(9) << "Si0";
	static vector<string> x(m+n+1); // Вектор Xi и F
	if (count == 0)
	{
		for (unsigned int i = {0}; i < m+n; ++i)
		{
			char temp[24]; 
			itoa(i+1,temp,10); // Переводим номер Икса в строку
			x[i] = (chk == 0) ? "X" : "Y";
			x[i] += temp;
		}
		x[m+n] = (chk == 0) ? "F " : "Ф ";
	}
	if (count > 0)
	{
		for (unsigned int i = {0}; i < n; ++i)
		{
			if (i == rj - 1)
			{
				for (unsigned int j = {0}; j < m+n; ++j)
				{
					if (j == ri)
					{
						string temp = x[i];
						x[i] = x[j+n];
						x[j+n] = temp; // Меняем X-ы, стоящие в разрешающей строке и столбце
					}
				}
			}
		}
	}
	count++;
	
	unsigned int k = {0};
	for (; k < n; ++k)
	{
		cout << std::setw(7) << x[k];
	}
	cout << endl;
	for (unsigned int i = {0}; i < m+1; ++i)
	{
		cout << x[k++];
		for (unsigned int j = {0}; j < n+1; ++j)
		{
			cout << std::setw(7) << std::fixed << std::setprecision(2) << arr[i][j];
		}
		cout << endl;
	}
	return x;
}

void set_new_data(std::vector<double> &c, std::vector<double> &b, double** a)
{
	using std::cout;
	using std::endl;

	chk++;
	count = 0;
	for (unsigned int i = {0}; i < c.size(); ++i)
			c[i] = -c[i];
	std::ofstream fout;
	fout.open("lab2.txt");
	if (!fout.is_open())
	{
		cout << endl << "Ошибка открытия файла.";
		exit(EXIT_FAILURE);
	}

	if (min_max == "max")
	{
		for (unsigned int i = {0}; i < c.size(); ++i)
			c[i] = -c[i];
		fout << "min" << endl;
	}
	else if (min_max == "min")
	{
		for (unsigned int i = {0}; i < c.size(); ++i)
			c[i] = -c[i];
		for (unsigned int i = {0}; i < b.size(); ++i) // Ф -> max --> -Ф = -by -> max
			b[i] = -b[i];
		fout << "max" << endl;
	}
	// Транспонируем матрицу системы ограничений A
	for (unsigned int i = {0}; i < b.size() / 2; ++i)
	{
		for (unsigned int j = {0}; j < c.size(); ++j)
		{
			if (i != j)
			{
				double temp = a[i][j];
				a[i][j] = a[j][i];
				a[j][i] = temp;
			}
		}
	}
	for (unsigned int i = b.size() - 1; i > b.size() / 2; --i)
	{
		for (unsigned int j = c.size() - 1; j > 0; --j)
		{
			if (i != j)
			{
				double temp = a[i][j];
				a[i][j] = a[j][i];
				a[j][i] = temp;
			}
		}
	}
	for (unsigned int i = {0}; i < b.size(); ++i)
	{
		for (unsigned int j = {0}; j < c.size(); ++j)
		{
			a[i][j] = -a[i][j];
		}
	}
	c.swap(b); // Меняем с и b, далее меняем c и b в файле lab2.txt
	// Перезаписываем значения для векторов c,b,a в файле - альтернативный метод транспонирования
	for (unsigned int i = {0}; i < c.size(); ++i)
		fout << c[i] << endl;
	fout << "q" << endl;
	for (unsigned int i = {0}; i < b.size(); ++i)
		fout << b[i] << endl;
	fout << "q" << endl;
	for (unsigned int i = {0}; i < b.size(); ++i)
	{
		for (unsigned int j = {0}; j < c.size(); ++j)
			fout << a[i][j] << endl;
	}
	fout.close();
	cout << endl << endl;
	show_new_task(c,b,a);
}

void show_new_task(std::vector<double> c, std::vector<double> b, double** a)
{
	using std::cout;
	using std::cin;
	using std::endl;

	cout << endl << "Постановка двойственной задачи (ДЗ) ЛП. Требуется найти решение следующей задачи:" 
	<< endl << std::setw(30) << "Ф = b^T * y -> " <<  ((min_max == "max") ? "min" : "max") << endl 
	<< std::setw(30) << "A^T * y >= c^T" << endl << std::setw(25) << "y >= 0" << endl
	<< "Каноническая форма ДЗ ЛП:" << endl;
	if (min_max == "max")
		cout << std::setw(17) << "Ф = ";
	else
		cout << std::setw(17) << "-Ф = ";
	for (unsigned int i = {0}; i < c.size(); ++i)
	{
		if (c[i] < 0)
			cout << c[i] << "y" << i+1;
		else
			cout << "+" << c[i] << "y" << i+1;
	}
	cout << " -> ";
	if (min_max == "min")
		cout <<  min_max << endl << endl; 
	else
		cout <<  "min" << endl << endl; 
	unsigned int k = c.size() + 1;
	for (unsigned int i = {0}; i < b.size(); ++i)
	{
		cout << "                 ";
		for (unsigned int j = {0}; j < c.size(); ++j)
		{
			if (a[i][j] != 0)
			{
				if (a[i][j] != 1 && a[i][j] != -1)
				{
					if (a[i][j] < 0 || j == 0)
						cout << a[i][j] << "y" << j+1;
					else
						if (a[i][j-1] != 0)
							cout << "+" << a[i][j] << "y" << j+1;
						else 
							cout << a[i][j] << "y" << j+1;
				}
				else
				{
					if (a[i][j] < 0)
						cout << "-y" << j+1;
					else if (j != 0)
						cout << "+" << "y" << j+1;
					else 
						cout << "y" << j+1;
				}
		
			}		
		}
		cout << "+y" << k++ << " = " << b[i] << endl;
	}
	cout << endl << "             ";
	for (unsigned int i = {0}; i < c.size()+b.size(); ++i)
	{
		if (i != c.size()+b.size()-1)
			cout << "y" << i+1 << ", ";
		else
			cout << "y" << i+1;
	}
	cout << " >= 0" << endl << "Алгоритм преобразования симплекс-таблицы (жордановы исключения):";

}