#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;

#define DX 0.25
#define DY 0.25

class Cell {
private:
	int index_x;
	int index_y;
	double value;
public:
	Cell() {};
	Cell(int x, int y, double v = 0.0) : index_x(x), index_y(y), value(v) {};
	int set_index_x(int ind_x) { index_x = ind_x; return index_x; };
	int set_index_y(int ind_y) { index_y = ind_y; return index_y; };
	double set_value(double v) { value = v; return value; };
	int get_index_x() { return index_x; };
	int get_index_y() { return index_y; };
	double get_value() { return value; };
};

class MeshField {
private:
	Cell **field;
	int width;
	int height;
public:
	MeshField(int n, int m);
	void showField();
	void initField();
	void setLoads(double ,double);
	int count_variable();
	Cell getCell(int i, int j) { return field[i][j]; };
	int get_width() { return width; };
	int get_height() { return height; };
	void initTempByResults(double*);
	int writeIntoFile(const char*);
};

class Variable {
private:
	int ind_x;
	int ind_y;
	double value;
	int pos;
public:
	Variable() { ind_x = ind_y = pos = 0; value = 0.0; };
	int get_pos() { return pos; };
	int get_ind_x() { return ind_x; };
	int get_ind_y() { return ind_y; };
	double get_value() { return value; };
	int set_pos(int p) { pos = p; return pos; };
	int set_ind_x(int x) { ind_x = x; return ind_x; };
	int set_ind_y(int y) { ind_y = y; return ind_y; };
	double set_value(double v) { value = v; return value; };
};

class Equation {
private:
	Variable **matx;
	int size;
public:
	Equation(int);
	int get_size() { return size; };
	void equationInit(MeshField);
	void showMatx();
	double *Gausse();
};

MeshField::MeshField(int n, int m) {
	width = m/DX + 2 + 1;
	height = n/DY + 2 + 1;
	field = new Cell*[height];
	for (int i = 0; i < height; i++) {
		field[i] = new Cell[width];
	}
}

void MeshField::initField() {
	for (int i = 0; i < height; i++) {
		for (int j = 0; j < width; j++) {
			field[i][j].set_index_x(i);
			field[i][j].set_index_y(j);
			field[i][j].set_value(0.0);
		}
	}

	for (int i = 1; i < height - 1; i++) {
		for (int j = i; j < width - i; j++)
			field[height - i - 1][j].set_value(-1.0);
	}

}

void MeshField::setLoads(double l, double b) {

	for (int j = 1; j < width - 1; j++)
		field[height - 2][j].set_value(b);

	for (int i = 1; i < height - 1; i++) {
		field[height - i - 1][i].set_value(l);
	}
}

int MeshField::count_variable() {
	int cv = 0;
	for (int i = 0; i < height; i++)
		for (int j = 0; j < width; j++)
			if (field[i][j].get_value() != 0.0)
				cv++;
	return cv;
}

void MeshField::showField() {
	cout << "..::: Field :::.." << endl;
	for (int i = 0; i < height; i++) {
		for (int j = 0; j < width; j++) {
			cout << field[i][j].get_value() << "\t";
		}
		cout << endl << endl;
	}

	cout << "__________________" << endl;
}

Equation::Equation(int n_var) {
	size = n_var;
	matx = new Variable*[n_var];
	for (int i = 0; i < n_var; i++)
		matx[i] = new Variable[n_var + 1];
}


void Equation::equationInit(MeshField mesh) {
	for (int i = 0; i < size; i++)
		for (int j = 0; j < size + 1; j++)
			matx[i][j].set_value(0.0);

	int pos = 0;
	/*init for every row*/
	for (int i = 0; i < size; i++) {
		for (int j = 1; j < mesh.get_height() - 1; j++) {
			for (int k = 1; k < mesh.get_width() - 1; k++) {
				if (mesh.getCell(j, k).get_value() != 0.0) {
					double v = mesh.getCell(j, k).get_value();
					matx[i][pos].set_value(v);
					matx[i][pos].set_pos(pos);
					matx[i][pos].set_ind_x(j);
					matx[i][pos].set_ind_y(k);
					pos++;
				}
			}
		}
		pos = 0;
	}
	/*end init for every row*/

	/*setting constants value*/
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			if (i == j) {
				if (matx[i][j].get_value() > 0.0) {
					matx[i][size].set_value(matx[i][size].get_value() + matx[i][j].get_value());
					matx[i][j].set_value(1.0);
					continue;
				}
			}
			matx[i][j].set_value(0.0);
		}
	}
	/*end setting const value*/


	/*setting boder 2nd type*/

	int lead_pos = -1;
	int sec_pos = -1;
	int equation_pos = -1;

	for (int i = 1; i < mesh.get_height() - 1; i++) {
		for (int j = 1; j < mesh.get_width() - 1; j++) {
			if (mesh.getCell(i, j).get_value() != 0.0)
				equation_pos++;
			if (mesh.getCell(i, j).get_value() == -1.0) {
				if (mesh.getCell(i, j + 1).get_value() == 0.0) {
					/*left side*/
					cout << "left" << endl;

					for (int k = 0; k < size; k++) {
						if (matx[equation_pos][k].get_ind_x() == i && matx[equation_pos][k].get_ind_y() == j) {
							//cout << "var pos = " << matx[equation_pos][k].get_pos() << endl;
							lead_pos = matx[equation_pos][k].get_pos();
						}
						if (matx[equation_pos][k].get_ind_x() == i + 1 && matx[equation_pos][k].get_ind_y() == j - 1) {
							//cout << "var scnd pos = " << matx[equation_pos][k].get_pos() << endl;
							sec_pos = matx[equation_pos][k].get_pos();
						}
						if (lead_pos != -1 && sec_pos != -1) {
							matx[equation_pos][lead_pos].set_value(1.0);
							matx[equation_pos][sec_pos].set_value(-1.0);
							matx[equation_pos][size].set_value(40.0 * sqrt(2.0)*DX);
							lead_pos = sec_pos = -1;
						}
					}

				}
				if (mesh.getCell(i, j - 1).get_value() == 0.0) {
					/*rigth side*/
					cout << "rigth" << endl;
				}
				if (mesh.getCell(i + 1, j).get_value() == 0.0) {
					/*bottom side*/
					cout << "bottom" << endl;
				}
			}
		}
	}
	equation_pos = -1;
	lead_pos = -1;
	/*end init 2nd type border*/

	/*init other dots*/
	int top_pos = -1;
	int rigth_pos = -1;
	int left_pos = -1;
	int bottom_pos = -1;

	for (int i = 1; i < mesh.get_height() - 1; i++) {
		for (int j = 1; j < mesh.get_width() - 1; j++) {
			if (mesh.getCell(i, j).get_value() != 0.0)
				equation_pos++;
			if (mesh.getCell(i, j).get_value() == -1.0 &&

				mesh.getCell(i, j - 1).get_value() != 0.0 &&
				mesh.getCell(i, j + 1).get_value() != 0.0 &&

				mesh.getCell(i - 1, j).get_value() != 0.0 &&
				mesh.getCell(i + 1, j).get_value() != 0.0) {
				
				for (int k = 0; k < size; k++) {
					if (matx[equation_pos][k].get_ind_x() == i && matx[equation_pos][k].get_ind_y() == j) {
						//cout << "hard equation" << endl;
						//cout << matx[equation_pos][k].get_pos() << endl;
						//cout << "eqtuation pos = " << equation_pos << endl;
						lead_pos = matx[equation_pos][k].get_pos();
					}
					if (matx[equation_pos][k].get_ind_x() == i && matx[equation_pos][k].get_ind_y() == j - 1) {
						//cout << matx[equation_pos][k].get_pos() << endl;
						left_pos = matx[equation_pos][k].get_pos();
					}
					if (matx[equation_pos][k].get_ind_x() == i && matx[equation_pos][k].get_ind_y() == j + 1) {
						//cout << matx[equation_pos][k].get_pos() << endl;
						rigth_pos = matx[equation_pos][k].get_pos();
					}
					if (matx[equation_pos][k].get_ind_x() == i - 1 && matx[equation_pos][k].get_ind_y() == j) {
						//cout << matx[equation_pos][k].get_pos() << endl;
						top_pos = matx[equation_pos][k].get_pos();
					}
					if (matx[equation_pos][k].get_ind_x() == i + 1 && matx[equation_pos][k].get_ind_y() == j) {
						//cout << matx[equation_pos][k].get_pos() << endl;
						bottom_pos = matx[equation_pos][k].get_pos();
					}
					if (lead_pos != -1 && left_pos != -1 && rigth_pos != -1 && top_pos != -1 && bottom_pos != -1) {
						matx[equation_pos][lead_pos].set_value(-4.0/DX);
						matx[equation_pos][left_pos].set_value(1.0/DX);
						matx[equation_pos][rigth_pos].set_value(1.0/DX);
						matx[equation_pos][top_pos].set_value(1.0/DY);
						matx[equation_pos][bottom_pos].set_value(1.0/DY);
						lead_pos = left_pos = rigth_pos = top_pos = bottom_pos = -1;
					}
				}
			}
		}
	}
	/*end init other dots*/
}

double *Equation::Gausse() {
	cout << endl;

	double **a;
	a = new double*[size];
	for (int i = 0; i < size; i++) {
		a[i] = new double[size];
	}

	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			a[i][j] = matx[i][j].get_value();
		}
	}

	double *y;

	y = new double[size];

	for (int i = 0; i < size; i++) {
		y[i] = matx[i][size].get_value();
	}

	int n = size;
	double *x, max;
	int k, index;
	const double eps = 0.0000000000000000000000000000000000000000000000000000000000000000000000000000001;  // точность
	x = new double[n];
	k = 0;
	while (k < n)
	{
		// Поиск строки с максимальным a[i][k]
		max = abs(a[k][k]);
		index = k;
		for (int i = k + 1; i < n; i++)
		{
			if (abs(a[i][k]) > max)
			{
				max = abs(a[i][k]);
				index = i;
			}
		}
		// Перестановка строк
		if (max < eps)
		{
			// нет ненулевых диагональных элементов
			cout << "Решение получить невозможно из-за нулевого столбца ";
			cout << index << " матрицы A" << endl;
			//return 0;
		}
		for (int j = 0; j < n; j++)
		{
			double temp = a[k][j];
			a[k][j] = a[index][j];
			a[index][j] = temp;
		}
		double temp = y[k];
		y[k] = y[index];
		y[index] = temp;
		// Нормализация уравнений
		for (int i = k; i < n; i++)
		{
			double temp = a[i][k];
			if (abs(temp) < eps) continue; // для нулевого коэффициента пропустить
			for (int j = 0; j < n; j++)
				a[i][j] = a[i][j] / temp;
			y[i] = y[i] / temp;
			if (i == k)  continue; // уравнение не вычитать само из себя
			for (int j = 0; j < n; j++)
				a[i][j] = a[i][j] - a[k][j];
			y[i] = y[i] - y[k];
		}
		k++;
	}
	// обратная подстановка
	for (k = n - 1; k >= 0; k--)
	{
		x[k] = y[k];
		for (int i = 0; i < k; i++)
			y[i] = y[i] - a[i][k] * x[k];
	}

	for (int i = 0; i < n; i++) {
		cout << x[i] << endl;
	}

	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			matx[i][j].set_value(abs(a[i][j]));
		}
	}
	for (int i = 0; i < size; i++){
		cout << matx[i][size].set_value(y[i]) << endl;
	}
	return y;
}

void MeshField::initTempByResults(double *y) {
	int count = 0;
	for (int i = 1; i < height - 1; i++) {
		for (int j = 1; j < width; j++) {
			if (field[i][j].get_value() != 0) {
				field[i][j].set_value(y[count]);
				count++;
			}
		}
	}
}

int MeshField::writeIntoFile(const char* fn){
	ofstream fout(fn);
	for(int i = 0; i < height; i++){
		for(int j = 0; j < width; j++){
			fout << i << "\t";
			fout << j << "\t";
			fout << field[i][j].get_value() << "\n";
		}
		fout << "\n";
	}
	return 1;
}

void Equation::showMatx() {

	for (int i = 0; i < size; i++) {
		cout << endl;
		for (int j = 0; j < size + 1; j++) {
			if (j == size) cout << "| ";
			cout << matx[i][j].get_value() << " "; 
		}
	}
}


int main() {
	MeshField meshFiled(4, 8);
	meshFiled.initField();
	meshFiled.setLoads(200.0, 100.0);
	meshFiled.showField();
	cout << meshFiled.count_variable();

	cout << endl;

	Equation equation(meshFiled.count_variable());
	equation.equationInit(meshFiled);
	equation.showMatx();

	cout << endl;
	double *y = equation.Gausse();

	meshFiled.initTempByResults(y);
	meshFiled.showField();
	meshFiled.setLoads(200.0, 100.0);
	meshFiled.writeIntoFile("temp.data");

	return 0;
}


