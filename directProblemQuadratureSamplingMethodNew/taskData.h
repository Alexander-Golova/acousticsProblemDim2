#pragma once

struct Point
{
	float x;
	float y;
};

// ������� ������������� ����
const float OMEGA = 1.0f;
const float C_0 = 1.0f;

// ���������� ����������
const float receiver = 11.0f;

// ���������� ����������� �� ������� ���������
const size_t NUMBER_PARTITION_POINTS = 50;

// ������ �������� � ������� ��������� ��������������
const float DOMAIN_IN_HOMOGENEITY = 10.0;
