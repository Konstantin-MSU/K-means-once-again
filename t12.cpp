// 212-Sushchev K-means

#include <iostream>     // библиотека вывода в консоль
#include <fstream>      // библиотека работы с файлами
#include <random>		// генератор сл. чисел
#include <vector>		// вектор
#include <algorithm>
#include <cmath>

class Interface {
	// 212-Sushchev - чтение файлов, введение параметров	
protected:

	std::ifstream V_file{ "variables.txt" }; // читаем файл с переменными
	std::vector <float> V; // вектор из всех параметров задачи  
	float w = 0; //переменная пробегающая все значения файла

	float Dens = 0; //число точек в 1 конусе
	float K = 0; //число конусов

	std::vector <float> X_vec; // координаты центров и радиусы
	std::vector <float> Y_vec;
	std::vector <float> R_vec;

	float Eps;

public:
	Interface() {
		if (!V_file.is_open()) {
			std::cout << "Error, file is not found" << std::endl;
			return;
		}
		while (V_file >> w) {
			V.push_back(w); // заполняем файл данными из листа
		}
		V_file.close();
		Dens = V[0];
		K = V[1];
		for (w = 0; w < K; w++) {
			X_vec.push_back(V[2 + w]);
			Y_vec.push_back(V[2 + w + K]);
			R_vec.push_back(V[2 + w + 2 * K]);
		}
		Eps = V[2 + 3 * K];
	}
};





class Tpair {
	// 212-Sushchev-класс работы с точкой, как парой координат со свойствами 
private:
	float X = 0;
	float Y = 0;
	int cell_X;
	int cell_Y;
	int Cl = 0;

public:
	Tpair() = default;
	Tpair(float x, float y) : X(x), Y(y) {}

	// сеттеры и геттеры
	void setX(float x) { X = x; }
	void setY(float y) { Y = y; }
	void set_cell_X(float x) { cell_X = x; }
	void set_cell_Y(float y) { cell_Y = y; }
	void setCl(int c) { Cl = c; }

	float getX() { return X; }
	float getY() { return Y; }
	int get_cell_X() { return cell_X; }
	int get_cell_Y() { return cell_Y; }
	int getCl() { return Cl; }

};



class Cone : public Interface {
	// 212-Sushchev-класс работы с конусом 
private:

	int Cone_Num; // порядковый номер конуса
	float coneX;
	float coneY;
	float coneR;
	int i = 0;

	std::vector <Tpair> inside;

public:
	Cone(int N) {
		Cone_Num = N; // присваеваем параметры в соответствии с номером
		coneX = X_vec[N];
		coneY = Y_vec[N];
		coneR = R_vec[N];

		for (i = 0; i < Dens; i++) {
			static std::random_device rd;
			static std::mt19937 seed(rd()); // создаем генератор с нормальным распределением
			std::normal_distribution<double> generated_X(coneX, coneR / 2);
			std::normal_distribution<double> generated_Y(coneY, coneR / 2);

			inside.push_back(Tpair(generated_X(seed), generated_Y(seed))); // заполняем маконус точками
		}
	}
	// геттеры
	float getCX() { return coneX; }
	float getCY() { return coneY; }
	float getCR() { return coneR; }
	std::vector<Tpair>& points() { return inside; }
};




class Pole : public Interface {
	//212-Sushchev класс самой плоскости на которой находятся точки
protected:
	std::vector<Tpair> Points; // множество точек
	float netR = V[3 + 3 * K]; // размер сетки (радиус) 
	int netA = V[4 + 3 * K]; // длина радиуса в клетках
	int K_Means = V[6 + 3 * K]; // K для K-means
	int i = 0;
	Tpair A = Tpair(0, 0);

public:
	Pole() {
		for (int i = 0; i < K; i++) { // перебираем точки из конусов, вставляем их в список всех точек
			Cone c(i);
			for (Tpair& A : c.points()) {

				A.set_cell_X(1 + int(A.getX() / netA));
				A.set_cell_Y(1 + int(A.getY() / netA));

				Points.push_back(A);

			}
		}
	}
	std::vector<Tpair>& getPoints() { return Points; }

};






class Distance : public Pole {
	// 212-Sushchev - класс работы с расстояниями между точками

public: Distance() {}

	  static float between(Tpair A, Tpair B) { // функция вычисления дистанции между 2 точками
		  return std::sqrt((A.getX() - B.getX()) * (A.getX() - B.getX()) +
			  (A.getY() - B.getY()) * (A.getY() - B.getY()));
	  }

	  Tpair closest(Tpair A) { // находит ближайщую точку
		  float D = 100000;
		  Tpair C = Tpair(0, 0);
		  for (Tpair B : Points) {
			  if ((between(A, B) < D) && (between(A, B) > Eps)) {
				  D = between(A, B);
				  C = B;
			  }
		  }
		  return C;
	  }
};







class Graph : public Distance {
	// 212-Sushchev - графф и вункции работы с ним
protected:
	float Max_edge = V[5 + 3 * K];
	int Edge_count = 0;
public:

	int Edge_Exist(Tpair A, Tpair B) { // проверка существования ребра
		if (between(A, B) < Eps) { return 0; }
		if (between(A, B) < Max_edge) { return 1; }
		return { 0 };
	}

	int Neigh(Tpair A) { // счет соседей у точки
		int D = 0;
		for (Tpair B : Points) {
			D += Edge_Exist(A, B);
		}
		return D;
	}

	int Dust(Tpair A) { // проверка, является ли точка аномалией (у нее меньше 3 соседей)
		if (Neigh(A) < 3) { return 0; }
		return 1;
	}

	Graph() {
		Edge_count = 0;
		for (Tpair A : Points) {
			Edge_count += Neigh(A);
		}
		Edge_count = int(Edge_count / 2);
	}
};






class cluster {
	// 212-Sushchev-класс работы с кластером
private:

	std::vector<Tpair> points; // точки кластера
	int points_num = 0;
	int cl_Num = 0; // нромер кластера
	Tpair Mean; // среднее значение

public:
	cluster() = default;
	void add(Tpair& A) { // добавление точки в кластер

		A.setCl(cl_Num);
		points.push_back(A);
		points_num += 1;
	}


	std::vector<Tpair> get_points() { return points; } // получение точек кластера

	void set_cl_Num(int a) { cl_Num = a; }
	int get_cl_Num() { return cl_Num; }

	Tpair get_Mean() { 
		Tpair B = Tpair(0, 0);
		for (Tpair A : points) {
			B = sum(A, B);
		}
		Mean = Tpair(B.getX() / points_num, B.getY() / points_num);
		
		return Mean; }

	Tpair sum(Tpair A, Tpair B) {
		return Tpair(A.getX() + B.getX(), A.getY() + B.getY());
	}
};




class K_Means_algorithm : public Graph {
	// 212-Sushchev - класс исполнения алгоритма K_means
protected:
	std::vector<Tpair> Set_Means; // центры тяжести
	std::vector<Tpair> Set_Means_Old; // центры тяжести предыдущей итерации
	std::vector<cluster> sectors; // кластеры
public:

	int belong(Tpair A) { // функция определяет, к какому уластеру относится точка
		int w = 0;
		float D = 100000;
		for (int i = 0; i < K_Means; i++) {
			if (between(Set_Means_Old[i], A) < D) {
				D = between(Set_Means_Old[i], A);
				w = i;
			}
		}
		return w;
	}

	K_Means_algorithm() {
		int j = 0;

		for (int i = 0; i < K_Means; i++) { // выбираем начальные точки
			static std::mt19937 gen(std::random_device{}());
			std::uniform_int_distribution<int> dist(0, K * Dens-1);
			//std::cout << "Points size: " << Points.size() << std::endl;
			Set_Means.push_back(Points[dist(gen)]);
		}
		
		while (true) {

			i = 0;
			j = 0;

			Set_Means_Old.clear();
			sectors.clear();
			for (int i = 0; i < K_Means; i++) { // сохраняем старые центры и создаем кластеры				
				Set_Means_Old.push_back(Set_Means[i]);		

				sectors.push_back(cluster());
			}

			Set_Means.clear();
			int w = 0;

			for (Tpair A : Points) { //распределяем точки
				w = belong(A);
				A.setCl(w);
				sectors[w].add(A);			
			}
			//i = 0;
			for (int i = 0; i < K_Means; i++) { //сравниваем среднее
				Set_Means.push_back(sectors[i].get_Mean());
							
				if (between(Set_Means[i], Set_Means_Old[i]) > Eps) {
					j = 1;
					//break;
				}
			}
			if (j == 0) { //если проверка была пройдена 5 раз, останавливаемся
				break;
			}
		}
	}
};





class plot_K_means : public K_Means_algorithm{
	// 212-Sushchev - класс вывода результатов в файл для GNU plot
private:
	std::ofstream output{ "clusters.txt" }; // открываем файл для записи
public:
	plot_K_means() {
		if (!output.is_open()) {
			std::cout << "Error, cannot open output file" << std::endl;
			return;
		}

		// перебираем все точки из в порядке возрастания кластера
		//int i = 0;
		for (int i = 0; i < K_Means; i++) {

			for (Tpair A : sectors[i].get_points()) {
				// выводим X, Y и номер кластера через пробел
				output << A.getX() << " " << A.getY() << " " << i+1 << "\n";
			}
		}
	}
};





int main() {
	//212-Sushchev - итоговый мэйн
	plot_K_means();
	return 0;
}