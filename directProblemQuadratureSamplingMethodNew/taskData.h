#pragma once

struct Point
{
	float x;
	float y;
};

// задание характеристик поля
const float OMEGA = 1.0f;
const float C_0 = 1.0f;

// координаты приемников
const float receiver = 11.0f;

// количество квадратиков по каждому измерению
const size_t NUMBER_PARTITION_POINTS = 50;

// размер квадрата в котором находится неоднородность
const float DOMAIN_IN_HOMOGENEITY = 10.0;
