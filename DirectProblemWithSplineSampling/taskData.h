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
const double RECEIVER = 1.1;

// количество квадратиков по каждому измерению
const size_t NUMBER_PARTITION_POINT = 20;
const size_t SPLITTING = 20;

// размер квадрата в котором находится неоднородность
const double DOMAIN_IN_HOMOGENEITY = 1.0;
const double h = DOMAIN_IN_HOMOGENEITY / NUMBER_PARTITION_POINT;

const double stepReceiver = 0.1 / NUMBER_PARTITION_POINT;
