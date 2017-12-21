#pragma once

struct Point
{
	float x;
	float y;
};

// задание характеристик поля
const float OMEGA = 6.2831853f;
const float C_0 = 1.0f;

// координаты приемников
const float receiver = 1.1f;

// количество квадратиков по каждому измерению
const size_t NUMBER_PARTITION_POINT = 40;
const size_t N_SQUARED = (NUMBER_PARTITION_POINT + 1) * (NUMBER_PARTITION_POINT + 1);

// размер квадрата в котором находится неоднородность
const float DOMAIN_IN_HOMOGENEITY = 1.0f;

// шаг по сетке
const float h = DOMAIN_IN_HOMOGENEITY / NUMBER_PARTITION_POINT;
