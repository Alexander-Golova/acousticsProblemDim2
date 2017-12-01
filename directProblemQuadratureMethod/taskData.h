#pragma once

struct Point
{
	double x;
	double y;
};

// задание характеристик поля
const double OMEGA = 1.0;
const double C_0 = 1.0;

// координаты приемников
const double receiver = 1.2;

// количество квадратиков по каждому измерению
const size_t NUMBER_PARTITION_POINT = 40;
const size_t N_SQUARED = (NUMBER_PARTITION_POINT + 1) * (NUMBER_PARTITION_POINT + 1);

// размер квадрата в котором находится неоднородность
const double DOMAIN_IN_HOMOGENEITY = 1.0;

// шаг по сетке
const double h = DOMAIN_IN_HOMOGENEITY / NUMBER_PARTITION_POINT;
