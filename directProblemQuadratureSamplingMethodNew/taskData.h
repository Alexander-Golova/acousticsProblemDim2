#pragma once

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

// количество квадратиков по каждому измерению
const size_t NUMBER_PARTITION_POINT = 50;

// размер квадрата в котором находится неоднородность
const double DOMAIN_IN_HOMOGENEITY = 10.0;
