#pragma once

const std::complex<double> I(0.0, 1.0);

struct Posize_t
{
	double x;
	double y;
};

// задание характеристик поля
const double OMEGA = 1.0;
const double C_0 = 1.0;

// координаты приемников
const double receiver = 11.0;

// количество источников
const size_t NUMBER_SOURCE = 5;

// количество квадратиков по каждому измерению
const size_t NUMBER_PARTITION_Posize = 50;

// размер квадрата в котором находится неоднородность
const double DOMAIN_IN_HOMOGENEITY = 10.0;

//функция Ханкеля
std::complex<double> Hankel(const double x);

// функция Грина
std::complex<double> G(const double x_1, const double x_2, const double y_1, const double y_2);
