//  Portions Copyright 2019
// Xuesong Zhou
//   If you help write or modify the code, please also list your names here.
//   The reason of having Copyright info here is to ensure all the modified version, as a whole, under the GPL 
//   and further prevent a violation of the GPL.

// More about "How to use GNU licenses for your own software"
// http://www.gnu.org/licenses/gpl-howto.html


extern void g_ProgramStop();

#include <iostream>
#include <fstream>
#include <list> 
#include <omp.h>
#include <algorithm>
#include <time.h>
#include <functional>
#include<stdio.h>   


#include <stack>
#include <string>
#include <vector>
#include <map>
#include <sstream>
using namespace std;
using std::string;
using std::ifstream;
using std::vector;
using std::map;
using std::istringstream;
template <typename T>

// some basic parameters setting

//Pls make sure the _MAX_K_PATH > Agentlite.cpp's g_number_of_K_paths+g_reassignment_number_of_K_paths and the _MAX_ZONE remain the same with .cpp's defination
#define _MAX_LABEL_COST 99999.0
#define _MAX_K_PATH 20 
#define _MAX_ZONE 70

#define _MAX_DEMANDTYPES 5 //because of the od demand store format,the MAX_demandtype must >=g_DEMANDTYPES.size()+1
#define _MAX_TIMEPERIODS 6

#define _MAX_TIMESLOT_PerPeriod 40

#define MIN_PER_TIMESLOT 15

// Linear congruential generator 
#define LCG_a 17364
#define LCG_c 0
#define LCG_M 65521  // it should be 2^32, but we use a small 16-bit number to save memory
#define _MAX_STRING_LINE 30000
#define FLAG_PEORIOD_TIMESLOT 0


FILE* g_pFileDebugLog = NULL;
FILE* g_pFileOutputLog = NULL;

extern FILE* g_pFileDebugLog;



//below shows where the functions used in Agentlite.cpp come from!
//Utility.cpp

#pragma warning(disable: 4244)  // stop warning: "conversion from 'int' to 'float', possible loss of data"


class CCSVParser
{
public:
	char Delimiter;
	bool IsFirstLineHeader;
	ifstream inFile;
	string mFileName;
	vector<string> LineFieldsValue;
	vector<string> Headers;
	map<string, int> FieldsIndices;

	vector<int> LineIntegerVector;

public:
	void  ConvertLineStringValueToIntegers()
	{
		LineIntegerVector.clear();
		for (unsigned i = 0; i < LineFieldsValue.size(); i++)
		{
			std::string si = LineFieldsValue[i];
			int value = atoi(si.c_str());

			if (value >= 1)
				LineIntegerVector.push_back(value);

		}
	}
	vector<string> GetHeaderVector()
	{
		return Headers;
	}

	int m_EmptyLineCount;
	bool m_bDataHubSingleCSVFile;
	string m_DataHubSectionName;
	bool m_bLastSectionRead;

	bool m_bSkipFirstLine;  // for DataHub CSV files

	CCSVParser(void)
	{
		Delimiter = ',';
		IsFirstLineHeader = true;
		m_bSkipFirstLine = false;
		m_bDataHubSingleCSVFile = false;
		m_bLastSectionRead = false;
		m_EmptyLineCount++;
	}

	~CCSVParser(void)
	{
		if (inFile.is_open()) inFile.close();
	}


	bool OpenCSVFile(string fileName, bool b_required)
	{
		mFileName = fileName;
		inFile.open(fileName.c_str());

		if (inFile.is_open())
		{
			if (IsFirstLineHeader)
			{
				string s;
				std::getline(inFile, s);
				vector<string> FieldNames = ParseLine(s);

				for (size_t i = 0;i < FieldNames.size();i++)
				{
					string tmp_str = FieldNames.at(i);
					size_t start = tmp_str.find_first_not_of(" ");

					string name;
					if (start == string::npos)
					{
						name = "";
					}
					else
					{
						name = tmp_str.substr(start);
						//			TRACE("%s,", name.c_str());
					}


					FieldsIndices[name] = (int)i;
				}
			}

			return true;
		}
		else
		{
			if (b_required)
			{

				cout << "File " << fileName << " does not exist. Please check." << endl;
				//g_ProgramStop();
			}
			return false;
		}
	}


	void CloseCSVFile(void)
	{
		inFile.close();
	}



	bool ReadRecord()
	{
		LineFieldsValue.clear();

		if (inFile.is_open())
		{
			string s;
			std::getline(inFile, s);
			if (s.length() > 0)
			{

				LineFieldsValue = ParseLine(s);

				return true;
			}
			else
			{

				return false;
			}
		}
		else
		{
			return false;
		}
	}

	vector<string> ParseLine(string line)
	{
		vector<string> SeperatedStrings;
		string subStr;

		if (line.length() == 0)
			return SeperatedStrings;

		istringstream ss(line);


		if (line.find_first_of('"') == string::npos)
		{

			while (std::getline(ss, subStr, Delimiter))
			{
				SeperatedStrings.push_back(subStr);
			}

			if (line.at(line.length() - 1) == ',')
			{
				SeperatedStrings.push_back("");
			}
		}
		else
		{
			while (line.length() > 0)
			{
				size_t n1 = line.find_first_of(',');
				size_t n2 = line.find_first_of('"');

				if (n1 == string::npos && n2 == string::npos) //last field without double quotes
				{
					subStr = line;
					SeperatedStrings.push_back(subStr);
					break;
				}

				if (n1 == string::npos && n2 != string::npos) //last field with double quotes
				{
					size_t n3 = line.find_first_of('"', n2 + 1); // second double quote

					//extract content from double quotes
					subStr = line.substr(n2 + 1, n3 - n2 - 1);
					SeperatedStrings.push_back(subStr);

					break;
				}

				if (n1 != string::npos && (n1 < n2 || n2 == string::npos))
				{
					subStr = line.substr(0, n1);
					SeperatedStrings.push_back(subStr);
					if (n1 < line.length() - 1)
					{
						line = line.substr(n1 + 1);
					}
					else // comma is the last char in the line string, push an empty string to the back of vector
					{
						SeperatedStrings.push_back("");
						break;
					}
				}

				if (n1 != string::npos && n2 != string::npos && n2 < n1)
				{
					size_t n3 = line.find_first_of('"', n2 + 1); // second double quote
					subStr = line.substr(n2 + 1, n3 - n2 - 1);
					SeperatedStrings.push_back(subStr);
					size_t idx = line.find_first_of(',', n3 + 1);

					if (idx != string::npos)
					{
						line = line.substr(idx + 1);
					}
					else
					{
						break;
					}
				}
			}

		}

		return SeperatedStrings;
	}

	template <class T> bool GetValueByFieldName(string field_name, T& value, bool NonnegativeFlag = true)
	{


		bool required_field = true;
		bool print_out = false;
		if (FieldsIndices.find(field_name) == FieldsIndices.end())
		{
			if (required_field)
			{
				cout << "Field " << field_name << " in file " << mFileName << " does not exist. Please check the file." << endl;

				g_ProgramStop();
			}
			return false;
		}
		else
		{
			if (LineFieldsValue.size() == 0)
			{
				return false;
			}

			int size = (int)(LineFieldsValue.size());
			if (FieldsIndices[field_name] >= size)
			{
				return false;
			}

			string str_value = LineFieldsValue[FieldsIndices[field_name]];

			if (str_value.length() <= 0)
			{
				return false;
			}

			istringstream ss(str_value);

			T converted_value;
			ss >> converted_value;

			if (/*!ss.eof() || */ ss.fail())
			{
				return false;
			}

			if (NonnegativeFlag && converted_value < 0)
				converted_value = 0;

			value = converted_value;
			return true;
		}
	}


	bool GetValueByFieldName(string field_name, string& value)
	{
		if (FieldsIndices.find(field_name) == FieldsIndices.end())
		{
			return false;
		}
		else
		{
			if (LineFieldsValue.size() == 0)
			{
				return false;
			}

			unsigned int index = FieldsIndices[field_name];
			if (index >= LineFieldsValue.size())
			{
				return false;
			}
			string str_value = LineFieldsValue[index];

			if (str_value.length() <= 0)
			{
				return false;
			}

			value = str_value;
			return true;
		}
	}

};



template <typename T>
T **AllocateDynamicArray(int nRows, int nCols)
{
	T **dynamicArray;

	dynamicArray = new (std::nothrow) T*[nRows];

	if (dynamicArray == NULL)
	{
		cout << "Error: insufficient memory.";
		g_ProgramStop();

	}

	for (int i = 0; i < nRows; i++)
	{
		dynamicArray[i] = new (std::nothrow) T[nCols];

		if (dynamicArray[i] == NULL)
		{
			cout << "Error: insufficient memory.";
			g_ProgramStop();
		}


	}

	return dynamicArray;
}

template <typename T>
void DeallocateDynamicArray(T** dArray, int nRows, int nCols)
{
	if (!dArray)
		return;

	for (int x = 0; x < nRows; x++)
	{
		delete[] dArray[x];
	}

	delete[] dArray;

}

template <typename T>
T ***Allocate3DDynamicArray(int nX, int nY, int nZ)
{
	T ***dynamicArray;

	dynamicArray = new (std::nothrow) T**[nX];

	if (dynamicArray == NULL)
	{
		cout << "Error: insufficient memory.";
		g_ProgramStop();
	}

	for (int x = 0; x < nX; x++)
	{
		if (x % 1000 == 0)
		{
			cout << "allocating 3D memory for " << x << endl;
		}


		dynamicArray[x] = new (std::nothrow) T*[nY];

		if (dynamicArray[x] == NULL)
		{
			cout << "Error: insufficient memory.";
			g_ProgramStop();
		}

		for (int y = 0; y < nY; y++)
		{
			dynamicArray[x][y] = new (std::nothrow) T[nZ];
			if (dynamicArray[x][y] == NULL)
			{
				cout << "Error: insufficient memory.";
				g_ProgramStop();
			}
		}
	}

	for (int x = 0; x < nX; x++)
		for (int y = 0; y < nY; y++)
			for (int z = 0; z < nZ; z++)
			{
				dynamicArray[x][y][z] = 0;
			}
	return dynamicArray;

}

template <typename T>
void Deallocate3DDynamicArray(T*** dArray, int nX, int nY)
{
	if (!dArray)
		return;
	for (int x = 0; x < nX; x++)
	{
		for (int y = 0; y < nY; y++)
		{
			delete[] dArray[x][y];
		}

		delete[] dArray[x];
	}

	delete[] dArray;

}



template <typename T>
T ****Allocate4DDynamicArray(int nM, int nX, int nY, int nZ)
{
	T ****dynamicArray;

	dynamicArray = new (std::nothrow) T***[nX];

	if (dynamicArray == NULL)
	{
		cout << "Error: insufficient memory.";
		g_ProgramStop();
	}
	for (int m = 0; m < nM; m++)
	{
		if (m % 100 == 0)
			cout << "allocating 4D memory for " << m << endl;

		dynamicArray[m] = new (std::nothrow) T**[nX];

		if (dynamicArray[m] == NULL)
		{
			cout << "Error: insufficient memory.";
			g_ProgramStop();
		}

		for (int x = 0; x < nX; x++)
		{
			dynamicArray[m][x] = new (std::nothrow) T*[nY];

			if (dynamicArray[m][x] == NULL)
			{
				cout << "Error: insufficient memory.";
				g_ProgramStop();
			}

			for (int y = 0; y < nY; y++)
			{
				dynamicArray[m][x][y] = new (std::nothrow) T[nZ];
				if (dynamicArray[m][x][y] == NULL)
				{
					cout << "Error: insufficient memory.";
					g_ProgramStop();
				}
			}
		}
	}
	return dynamicArray;

}

template <typename T>
void Deallocate4DDynamicArray(T**** dArray, int nM, int nX, int nY)
{
	if (!dArray)
		return;
	for (int m = 0; m < nM; m++)
	{
		for (int x = 0; x < nX; x++)
		{
			for (int y = 0; y < nY; y++)
			{
				delete[] dArray[m][x][y];
			}

			delete[] dArray[m][x];
		}
		delete[] dArray[m];
	}
	delete[] dArray;

}


//struct MyException : public exception {
//	const char * what() const throw () {
//		return "C++ Exception";
//	}
//};
//

class CDemand_Period {
public:

	CDemand_Period()
	{
		demand_period_id = 0;
		starting_time_slot_no = 0;
		ending_time_slot_no = 0;

	}
	string demand_period;
	string time_period;
	int demand_period_id;
	int starting_time_slot_no;
	int ending_time_slot_no;

};


class CDemand_Type {
public:
	CDemand_Type()
	{
		demand_type_id = 0;
		value_of_time = 10;
		capacity_equivalent_ratio = 1.0;
	}

	string demand_type;
	int demand_type_id;
	float value_of_time;  // dollar per hour
	float capacity_equivalent_ratio;
};

class CDemand_Info {
public:
	float volume;
	float cost;
	float time;
	float distance;
	vector<int> path_node_vector;
	CDemand_Info()
	{
		volume = 0;
		cost = 0;
		time = 0;
		distance = 0;
	}
};

class Assignment {
public:
	Assignment()
	{
		g_demand_array = NULL;
		g_origin_demand_array = NULL;
		//pls check following 7 settings before running programmer
		g_number_of_threads = 4000;
		g_number_of_K_paths = 2;
		g_number_of_demand_periods = 24;
		g_reassignment_tau0 = 999;

		g_number_of_links = 0;
		g_number_of_nodes = 0;
		g_number_of_zones = 0;
		g_number_of_demand_types = 0;

		b_debug_detail_flag = 1;

	}

	void InitializeDemandMatrix(int number_of_zones, int number_of_demand_types)
	{
		g_number_of_zones = number_of_zones;
		g_number_of_demand_types = number_of_demand_types;

		g_demand_array = Allocate4DDynamicArray<CDemand_Info>(number_of_zones, number_of_zones, number_of_demand_types, _MAX_TIMEPERIODS);
		g_origin_demand_array = Allocate3DDynamicArray<float>(number_of_zones, number_of_demand_types, _MAX_TIMEPERIODS);


		for (int i = 0;i < number_of_zones;i++)
		{
			for (int dt = 0;dt < number_of_demand_types;dt++)
			{
				for (int tau = 0;tau < g_number_of_demand_periods;tau++)
				{

					g_origin_demand_array[i][dt][tau] = 0.0;
				}
			}

		}
		total_demand_volume = 0.0;
		for (int i = 0;i < number_of_demand_types;i++)
		{
			for (int tau = 0;tau < g_number_of_demand_periods;tau++)
			{
				total_demand[i][tau] = 0.0;
			}
		}

		g_DemandGlobalMultiplier = 1.0f;

	};
	~Assignment()
	{

		if (g_demand_array != NULL)
			Deallocate4DDynamicArray(g_demand_array, g_number_of_zones, g_number_of_zones, g_number_of_demand_types);
		if (g_origin_demand_array != NULL)
			Deallocate3DDynamicArray(g_origin_demand_array, g_number_of_zones, g_number_of_demand_types);
	}
	int g_number_of_threads;
	int g_number_of_K_paths;

	int g_reassignment_tau0;

	int b_debug_detail_flag;
	std::map<string, int> g_internal_node_to_seq_no_map;  // hush table, map external node number to internal node sequence no. 
	std::map<int, int> g_zoneid_to_zone_seq_no_mapping;// from integer to integer map zone_id to zone_seq_no

	CDemand_Info**** g_demand_array;
	float*** g_origin_demand_array;

	//StatisticOutput.cpp
	float total_demand_volume;
	//NetworkReadInput.cpp and ShortestPath.cpp



	std::vector<CDemand_Period> g_DemandPeriodVector;
	std::vector<CDemand_Type> g_DemandTypeVector;

	std::map<string, int> demand_period_to_seqno_mapping;
	std::map<string, int> demand_type_2_seqno_mapping;


	float total_demand[_MAX_DEMANDTYPES][_MAX_TIMEPERIODS];
	float g_DemandGlobalMultiplier;

	int g_number_of_links;
	int g_number_of_nodes;
	int g_number_of_zones;
	int g_number_of_demand_types;
	int g_number_of_demand_periods;

};

Assignment assignment;

class CVDF_Period
{
public:

	CVDF_Period()
	{
		m = 0.5;
	}


	int starting_time_slot_no;  // in 15 min slot
	int ending_time_slot_no;
	string period;


	//standard BPR parameter 
	float alpha;
	float beta;
	float capacity;
	float FFTT;

	//updated BPR-X parameters
	float gamma;
	float mu;
	float m;

	// inpput
	float volume;

	//output
	float avg_delay;
	float avg_travel_time = 0;
	float avg_waiting_time = 0;

	//float Q[_MAX_TIMESLOT_PerPeriod];  // t starting from starting_time_slot_no if we map back to the 24 hour horizon 
	float waiting_time[_MAX_TIMESLOT_PerPeriod];
	float arrival_rate[_MAX_TIMESLOT_PerPeriod];

	float get_waiting_time(int relative_time_slot_no)
	{
		if (relative_time_slot_no >=0 && relative_time_slot_no < _MAX_TIMESLOT_PerPeriod)
			return waiting_time[relative_time_slot_no];
		else
			return 0;

	}
	int t0, t3;

	void Setup()
	{

	}

	float  PeroformBPR(float volume)
	{

		return FFTT * (1 + alpha * pow(volume / max(0.00001, capacity), beta));

		// volume --> avg_traveltime

	}

	float PeroformBPR_X(float volume)
	{
		// Step 1: Initialization
		int L = ending_time_slot_no - starting_time_slot_no;  // in 15 min slot

		if (L >= _MAX_TIMESLOT_PerPeriod - 1)
			return 0;

		float mid_time_slot_no = starting_time_slot_no + L / 2.0;  // t1;
		for (int t = 0; t <= L; t++)
		{
			waiting_time[t] = 0;
			arrival_rate[t] = 0;
		}
		//int L = ending_time_slot_no - starting_time_slot_no;  // in 15 min slot

		// Case 1
		if (volume <= L * mu / 2)
		{
			float P = 0;
			t0 = mid_time_slot_no - P / 2.0;
			t3 = mid_time_slot_no + P / 2.0;
			int t2 = m * (t3 - t0) + t0;
			for (int t = 0; t <= L; t++)
			{
				waiting_time[t] = 0;
				arrival_rate[t] = (volume / L);
			}
			avg_waiting_time = 0;
			//cout << avg_waiting_time << endl;
			avg_travel_time = FFTT + avg_waiting_time;
		}

		// Case 2
		if (volume > L * mu / 2 && volume <= (L * mu / 2) + L)
		{
			float P = volume * 2 / mu - L;
			t0 = mid_time_slot_no - P / 2.0;
			t3 = mid_time_slot_no + P / 2.0;
			int t2 = m * (t3 - t0) + t0;
			for (int tt = 0; tt <= L; tt++)
			{
				int time = starting_time_slot_no + tt;
				if (time < t0)
				{
					waiting_time[tt] = 0;
					arrival_rate[tt] = mu / 2;
				}
				if (time > t0 && time <= t3)
				{
					waiting_time[tt] = 1 / (4.0*mu) *gamma *(time - t0)*(time - t0) * (time - t3)*(time - t3);
					arrival_rate[tt] = gamma * (time - t0)*(time - t2)*(time - t3) + mu;
				}
				if (time > t3)
				{
					waiting_time[tt] = 0;
					arrival_rate[tt] = mu / 2;
				}
				avg_waiting_time = gamma / (120 * mu)*pow(P, 4.0);
				//cout << avg_waiting_time << endl;
				avg_travel_time = FFTT + avg_waiting_time;
			}
		}

		// Case 3
		if (volume > (L * mu / 2) + L)
		{
			float P = volume * 2 / mu - L;
			t0 = mid_time_slot_no - P / 2.0;
			t3 = mid_time_slot_no + P / 2.0;
			int t2 = m * (t3 - t0) + t0;
			for (int tt = 0; tt <= L; tt++)
			{
				int time = starting_time_slot_no + tt;
				waiting_time[tt] = 1 / (4.0*mu) *gamma *(time - t0)*(time - t0) * (time - t3)*(time - t3);
				arrival_rate[tt] = gamma * (time - t0)*(time - t2)*(time - t3) + mu;
				if (arrival_rate[tt] < 0)
					avg_waiting_time = gamma / (120 * mu)*pow(P, 4.0);
				//cout << avg_waiting_time << endl;
				avg_travel_time = FFTT + avg_waiting_time;
			}
		}
		//cout << avg_travel_time << endl;
		return avg_travel_time;

		//float P1 = max(0,volume * 2 / mu - L);
		//float P = min(P1, L);
		////float mid_time_slot_no = starting_time_slot_no + L / 2.0;  // t1;
		//t0 = mid_time_slot_no - P / 2.0;
		//t3 = mid_time_slot_no + P / 2.0;

		//for (int t = 0; t <= L; t++)
		//{
		//	waiting_time[t] = 0;
		//	arrival_rate[t] = 0;
		//}

		//int t;

		//for (int tt = 0; tt <= P; tt++)
		//{
		//	t = t0 + tt;
		//	waiting_time[t- starting_time_slot_no] = 1 / (4.0*mu) *gamma *(t - t0)*(t - t0) * (t - t3)*(t - t3);
		//}

		//if (volume>0)
		//{ 
		//	int t2 = m * (t3 - t0) + t0;
		//	for (int tt = 0; tt <= L; tt++)
		//	{
		//		t = tt + starting_time_slot_no;
		//		arrival_rate[t - starting_time_slot_no] = gamma * (t - t0)*(t - t2)*(t - t3) + mu;
		//	}
		//}
		//

		//float avg_w = gamma / (120 * mu)*pow(P, 4.0) * (P / L) + 0 * (L - P) / L;
		//return (FFTT + avg_w);
		//
	}

};


class CLink
{
public:
	CLink()  // construction 
	{

		free_flow_travel_time_in_min = 1;

		for (int tau = 0; tau < _MAX_TIMEPERIODS; tau++)
		{
			flow_volume_per_period[tau] = 0;
			queue_length_perslot[tau] = 0;
			travel_time_per_period[tau] = 0;


			TDBaseTT[tau] = 0;
			TDBaseCap[tau] = 0;
			TDBaseFlow[tau] = 0;
			TDBaseQueue[tau] = 0;


			//cost_perhour[tau] = 0;
		}
		link_spatial_capacity = 100;


	}

	~CLink()
	{
		//if (flow_volume_for_each_o != NULL)
		//	delete flow_volume_for_each_o;
	}

	void free_memory()
	{
	}

	void AddAgentsToLinkVolume()
	{


	}

	std::vector<int> m_total_volume_vector;
	std::vector<float> m_avg_travel_time;

	// 1. based on BPR. 


	int m_LeftTurn_link_seq_no;

	int m_RandomSeed;
	int link_seq_no;
	string link_id;
	int from_node_seq_no;
	int to_node_seq_no;


	float fftt;
	float free_flow_travel_time_in_min;


	CVDF_Period VDF_period[_MAX_TIMEPERIODS];

	float TDBaseTT[_MAX_TIMEPERIODS];
	float TDBaseCap[_MAX_TIMEPERIODS];
	float TDBaseFlow[_MAX_TIMEPERIODS];
	float TDBaseQueue[_MAX_TIMEPERIODS];


	int type;
	float link_spatial_capacity;

	//static
	//float flow_volume;
	//float travel_time;

	float flow_volume_per_period[_MAX_TIMEPERIODS];
	float queue_length_perslot[_MAX_TIMEPERIODS];  // # of vehicles in the vertical point queue
	float travel_time_per_period[_MAX_TIMEPERIODS];


	string demand_type_code;

	int number_of_periods;

	float length;
	//std::vector <SLinkMOE> m_LinkMOEAry;
	//beginning of simulation data 

	//toll related link
	//int m_TollSize;
	//Toll *pTollVector;  // not using SLT here to avoid issues with OpenMP
	std::map<int, float> TollMAP;

	void CalculateTD_VDFunction();

	float get_VOC_ratio(int tau)
	{

		return (flow_volume_per_period[tau] + TDBaseFlow[tau]) / max(0.00001, TDBaseCap[tau]);
	}

	float get_speed(int tau)
	{
		return length / max(travel_time_per_period[tau], 0.0001) * 60;  // per hour
	}


};


class CNode
{
public:
	CNode()
	{
		zone_id = -1;
		//accessible_node_count = 0;
	}

	//int accessible_node_count;

	int node_seq_no;  // sequence number 
	string node_id;      //external node number 
	int zone_id = -1;
	double x;
	double y;

	std::vector<CLink> m_outgoing_link_vector;

};


extern std::vector<CNode> g_node_vector;
extern std::vector<CLink> g_link_vector;


class COZone
{
public:
	int zone_seq_no;  // 0, 1, 
	int zone_id;  // external zone id // this is origin zone
	int node_seq_no;


};

extern std::vector<COZone> g_zone_vector;
extern std::map<int, int> g_zoneid_to_zone_seq_no_mapping;

class CAGBMAgent
{
public:

	int agent_id;
	int income;
	int gender;
	int vehicle;
	int purpose;
	int flexibility;
	float preferred_arrival_time;
	float travel_time_in_min;
	float free_flow_travel_time;
	int from_zone_seq_no;
	int to_zone_seq_no;
	int type;
	int time_period;
	int k_path;
	float volume;
	float arrival_time_in_min;



};
extern std::vector<CAGBMAgent> g_agbmagent_vector;


class NetworkForSP  // mainly for shortest path calculation
{
public:

	NetworkForSP()
	{
		//pFileAgentPathLog = NULL;
	}

	int m_threadNo;  // internal thread number 

	int m_ListFront; // used in coding SEL
	int m_ListTail;  // used in coding SEL
	int* m_SENodeList; // used in coding SEL

	float* m_node_label_cost;  // label cost // for shortest path calcuating
	float* m_label_time_array;  // time-based cost
	float* m_label_distance_array;  // distance-based cost

	int * m_node_predecessor;  // predecessor for nodes
	int* m_node_status_array; // update status 
	int* m_link_predecessor;  // predecessor for this node points to the previous link that updates its label cost (as part of optimality condition) (for easy referencing)

	int  m_iteration_k;
	int  m_origin_node; // assigned nodes for computing 
	int  m_origin_zone_seq_no;
	int  m_demand_time_period_no; // assigned nodes for computing 
	int  m_demand_type_no; // assigned nodes for computing 

	// major function 1:  allocate memory and initialize the data 
	void AllocateMemory(int number_of_nodes)
	{

		m_SENodeList = new int[number_of_nodes];  //1
		m_node_status_array = new int[number_of_nodes];  //2
		m_label_time_array = new float[number_of_nodes];  //3
		m_label_distance_array = new float[number_of_nodes];  //4
		m_node_predecessor = new int[number_of_nodes];  //5
		m_link_predecessor = new int[number_of_nodes];  //6
		m_node_label_cost = new float[number_of_nodes];  //7
	}

	~NetworkForSP()
	{
		if (m_SENodeList != NULL)  //1
			delete m_SENodeList;

		if (m_node_status_array != NULL)  //2
			delete m_node_status_array;

		if (m_label_time_array != NULL)  //3
			delete m_label_time_array;

		if (m_label_distance_array != NULL)  //4
			delete m_label_distance_array;

		if (m_node_predecessor != NULL)  //5
			delete m_node_label_cost;

		if (m_link_predecessor != NULL)  //6
			delete m_link_predecessor;

		if (m_node_label_cost != NULL)  //7
			delete m_node_label_cost;




		if (m_link_predecessor != NULL)
			delete m_link_predecessor;

	}


	// SEList: scan eligible List implementation: the reason for not using STL-like template is to avoid overhead associated pointer allocation/deallocation
	void SEList_clear()
	{
		m_ListFront = -1;
		m_ListTail = -1;
	}

	void SEList_push_front(int node)
	{
		if (m_ListFront == -1)  // start from empty
		{
			m_SENodeList[node] = -1;
			m_ListFront = node;
			m_ListTail = node;
		}
		else
		{
			m_SENodeList[node] = m_ListFront;
			m_ListFront = node;
		}
	}
	void SEList_push_back(int node)
	{
		if (m_ListFront == -1)  // start from empty
		{
			m_ListFront = node;
			m_ListTail = node;
			m_SENodeList[node] = -1;
		}
		else
		{
			m_SENodeList[m_ListTail] = node;
			m_SENodeList[node] = -1;
			m_ListTail = node;
		}
	}

	bool SEList_empty()
	{
		return(m_ListFront == -1);
	}

	int SEList_front()
	{
		return m_ListFront;
	}

	void SEList_pop_front()
	{
		int tempFront = m_ListFront;
		m_ListFront = m_SENodeList[m_ListFront];
		m_SENodeList[tempFront] = -1;
	}


	//major function: update the cost for each node at each SP tree, using a stack from the origin structure 

	int calculate_TD_link_flow(Assignment& assignment, int iteration_number);

	//major function 2: // time-dependent label correcting algorithm with double queue implementation
	int optimal_label_correcting(Assignment& assignment, int iteration_k)

	{
		int shortest_path_debugging_flag = 0;
		FILE* pFileDebugLog = NULL;
		if (m_iteration_k != iteration_k)
			return 0;

		int origin_node = m_origin_node; // assigned nodes for computing 
		int demand_type = m_demand_type_no; // assigned nodes for computing 

		int internal_debug_flag = 0;
		if (g_node_vector[origin_node].m_outgoing_link_vector.size() == 0)
		{
			return 0;
		}

		for (int i = 0; i < assignment.g_number_of_nodes; i++) //Initialization for all non-origin nodes
		{
			m_node_status_array[i] = 0;  // not scanned
			m_node_label_cost[i] = _MAX_LABEL_COST;
			//m_node_label_cost_withouttoll[i] = _MAX_LABEL_COST;
			m_label_time_array[i] = 0;
			m_label_distance_array[i] = 0;
			m_node_predecessor[i] = -1;  // pointer to previous NODE INDEX from the current label at current node and time
			m_link_predecessor[i] = -1;  // pointer to previous NODE INDEX from the current label at current node and time
		}

		//Initialization for origin node at the preferred departure time, at departure time, cost = 0, otherwise, the delay at origin node
		m_label_time_array[origin_node] = 0;
		m_node_label_cost[origin_node] = 0.0;
		//m_node_label_cost_withouttoll[origin_node] = 0.0;

		SEList_clear();
		SEList_push_back(origin_node);

		clock_t start_t = clock();
		while (!SEList_empty())
		{
			int from_node = SEList_front();//pop a node FromID for scanning

			SEList_pop_front();  // remove current node FromID from the SE list
			m_node_status_array[from_node] = 2;

			//if (shortest_path_debugging_flag)
			//	fprintf(pFileDebugLog, "SP: SE node: %d\n", g_node_vector[from_node].node_id);

			//scan all outbound nodes of the current node
			for (int i = 0; i < g_node_vector[from_node].m_outgoing_link_vector.size(); i++)  // for each link (i,j) belong A(i)
			{

				int to_node = g_node_vector[from_node].m_outgoing_link_vector[i].to_node_seq_no;

				/*if (to_node == origin_node)
					continue;*/

				bool  b_node_updated = false;
				/*if (g_link_vector[g_node_vector[from_node].m_outgoing_link_vector[i].link_seq_no].TollMAP[assignment.g_DemandTypeVector[demand_type]] / max(0.0001, assignment.g_VOT_PerDemandType_MAP[assignment.g_DemandTypeVector[demand_type]]) * 60.0f > 0)
				{
					int a = g_link_vector[g_node_vector[from_node].m_outgoing_link_vector[i].link_seq_no].TollMAP[assignment.g_DemandTypeVector[demand_type]] / max(0.0001, assignment.g_VOT_PerDemandType_MAP[assignment.g_DemandTypeVector[demand_type]]) * 60.0f;
					cout << "a" << endl;
				}*/

				float new_time = m_label_time_array[from_node] + g_link_vector[g_node_vector[from_node].m_outgoing_link_vector[i].link_seq_no].travel_time_per_period[m_demand_time_period_no];
				float new_distance = m_label_distance_array[from_node] + g_link_vector[g_node_vector[from_node].m_outgoing_link_vector[i].link_seq_no].length;

				float new_to_node_cost = m_node_label_cost[from_node] + g_link_vector[g_node_vector[from_node].m_outgoing_link_vector[i].link_seq_no].travel_time_per_period[m_demand_time_period_no] + 0;

				//g_link_vector[g_node_vector[from_node].m_outgoing_link_vector[i].link_seq_no].TollMAP[demand_type]] / max(0.0001, assignment.g_DemandTypeVector [demand_type]);
								//float new_to_node_cost_withouttoll = m_node_label_cost_withouttoll[from_node] + m_link_cost_withouttoll_array[link_entering_time_interval][demand_type][g_node_vector[from_node].m_outgoing_link_vector[i].link_seq_no];
								/*if (shortest_path_debugging_flag)
								{
									fprintf(pFileDebugLog, "SP: checking from node %d, to node %d  cost = %d\n",
										g_node_vector[from_node].node_id,
										g_node_vector[to_node].node_id,
										new_to_node_cost, g_node_vector[from_node].m_outgoing_link_vector[i].cost);
								}*/

				if (new_to_node_cost < m_node_label_cost[to_node]) // we only compare cost at the downstream node ToID at the new arrival time t
				{

					if (shortest_path_debugging_flag)
					{
						fprintf(pFileDebugLog, "SP: updating node: %d current cost: %.2f, new cost %.2f\n",
							g_node_vector[to_node].node_id,
							m_node_label_cost[to_node], new_to_node_cost);
					}

					// update cost label and node/time predecessor
					m_label_time_array[to_node] = new_time;
					m_label_distance_array[to_node] = new_distance;

					m_node_label_cost[to_node] = new_to_node_cost;
					//m_node_label_cost_withouttoll[to_node] = new_to_node_cost_withouttoll;
					int link_seq_no = g_node_vector[from_node].m_outgoing_link_vector[i].link_seq_no;

					m_node_predecessor[to_node] = from_node;  // pointer to previous physical NODE INDEX from the current label at current node and time
					m_link_predecessor[to_node] = g_node_vector[from_node].m_outgoing_link_vector[i].link_seq_no;  // pointer to previous physical NODE INDEX from the current label at current node and time

					b_node_updated = true;

					//if (shortest_path_debugging_flag)
					//	fprintf(pFileDebugLog, "SP: add node %d into SE List\n",
					//		g_node_vector[to_node].node_id);

					//to_node is zone centroid and not origin_node,is to make sure no passing zones, only needed in network with connector
					/*if (g_node_vector[to_node].zone_id != -1 && to_node != origin_node)
					{
						m_node_status_array[to_node] = 1;
					}*/

					if (m_node_status_array[to_node] == 0)
					{
						SEList_push_back(to_node);
						m_node_status_array[to_node] = 1;
					}
					if (m_node_status_array[to_node] == 2)
					{
						SEList_push_front(to_node);
						m_node_status_array[to_node] = 1;
					}

				}

			}
		}

		clock_t end_t = clock();
		clock_t total_t = (end_t - start_t);
		//		cout << "( total time of shortest path calculation = " << total_t << " milliseconds" << " " << ")" << endl;



		return 1;  // one to all shortest pat
	}


};

std::vector<CNode> g_node_vector;
std::vector<CLink> g_link_vector;
std::vector<COZone> g_zone_vector;
std::vector<CAGBMAgent> g_agbmagent_vector;
std::vector<NetworkForSP*> g_NetworkForSP_vector;

//    This file is part of FlashDTA.

//    FlashDTA is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.

//    FlashDTA  is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.

//    You should have received a copy of the GNU General Public License
//    along with DTALite.  If not, see <http://www.gnu.org/licenses/>.


int g_read_integer(FILE* f, bool speicial_char_handling)
// read an integer from the current pointer of the file, skip all spaces
{
	char ch, buf[32];
	int i = 0;
	int flag = 1;
	/* returns -1 if end of file is reached */

	while (true)
	{
		ch = getc(f);
		//cout << "get from node successful: " << ch;
		if (ch == EOF || (speicial_char_handling && (ch == '*' || ch == '$')))
			return -1; // * and $ are special characters for comments
		if (isdigit(ch))
			break;
		if (ch == '-')
			flag = -1;
		else
			flag = 1;
	};
	if (ch == EOF) return -1;


	while (isdigit(ch)) {
		buf[i++] = ch;
		//cout << "isdigit" << buf[i++] << endl;
		ch = fgetc(f);
		//cout << "new ch" << ch;
	}
	buf[i] = 0;


	return atoi(buf) * flag;

}


float g_read_float(FILE *f)
//read a floating point number from the current pointer of the file,
//skip all spaces

{
	char ch, buf[32];
	int i = 0;
	int flag = 1;

	/* returns -1 if end of file is reached */

	while (true)
	{
		ch = getc(f);
		if (ch == EOF || ch == '*' || ch == '$') return -1;
		if (isdigit(ch))
			break;

		if (ch == '-')
			flag = -1;
		else
			flag = 1;

	};
	if (ch == EOF) return -1;
	while (isdigit(ch) || ch == '.') {
		buf[i++] = ch;
		ch = fgetc(f);

	}
	buf[i] = 0;

	/* atof function converts a character string (char *) into a doubleing
	pointer equivalent, and if the string is not a floting point number,
	a zero will be return.
	*/

	return (float)(atof(buf) * flag);

}



//split the string by "_"
vector<string> split(const string &s, const string &seperator) {
	vector<string> result;
	typedef string::size_type string_size;
	string_size i = 0;

	while (i != s.size()) {
		int flag = 0;
		while (i != s.size() && flag == 0) {
			flag = 1;
			for (string_size x = 0; x < seperator.size(); ++x)
				if (s[i] == seperator[x]) {
					++i;
					flag = 0;
					break;
				}
		}

		flag = 0;
		string_size j = i;
		while (j != s.size() && flag == 0) {
			for (string_size x = 0; x < seperator.size(); ++x)
				if (s[j] == seperator[x]) {
					flag = 1;
					break;
				}
			if (flag == 0)
				++j;
		}
		if (i != j) {
			result.push_back(s.substr(i, j - i));
			i = j;
		}
	}
	return result;
}

vector<float> g_time_parser(vector<string>& inputstring)
{
	vector<float> output_global_minute;

	for (int k = 0; k < inputstring.size(); k++)
	{
		vector<string> sub_string = split(inputstring[k], "_");

		for (int i = 0; i < sub_string.size(); i++)
		{
			//HHMM
			//012345
			char hh1 = sub_string[i].at(0);
			char hh2 = sub_string[i].at(1);
			char mm1 = sub_string[i].at(2);
			char mm2 = sub_string[i].at(3);

			float hhf1 = ((float)hh1 - 48);
			float hhf2 = ((float)hh2 - 48);
			float mmf1 = ((float)mm1 - 48);
			float mmf2 = ((float)mm2 - 48);

			float hh = hhf1 * 10 * 60 + hhf2 * 60;
			float mm = mmf1 * 10 + mmf2;
			float global_mm_temp = hh + mm;
			output_global_minute.push_back(global_mm_temp);
		}
	}

	return output_global_minute;
} // transform hhmm to minutes 



void g_ProgramStop()
{

	cout << "AgentLite Program stops. Press any key to terminate. Thanks!" << endl;
	getchar();
	exit(0);
};



//void ReadLinkTollScenarioFile(Assignment& assignment)
//{
//
//	for (unsigned li = 0; li < g_link_vector.size(); li++)
//	{
//
//		g_link_vector[li].TollMAP.erase(g_link_vector[li].TollMAP.begin(), g_link_vector[li].TollMAP.end()); // remove all previouly read records
//	}
//
//	// generate toll based on demand type code in input_link.csv file
//	int demand_mode_type_count = 0;
//
//	for (unsigned li = 0; li < g_link_vector.size(); li++)
//	{
//		if (g_link_vector[li].demand_type_code.size() >= 1)
//		{  // with data string
//
//			std::string demand_type_code = g_link_vector[li].demand_type_code;
//
//			vector<float> TollRate;
//			for (int dt = 0; dt < assignment.g_DemandTypeVector.size(); dt++)
//			{
//				CString number;
//				number.Format(_T("%d"), dt);
//
//				std::string str_number = CString2StdString(number);
//				if (demand_type_code.find(str_number) == std::string::npos)   // do not find this number
//				{
//					g_link_vector[li].TollMAP[dt] = 999;
//					demand_mode_type_count++;
//				}
//				else
//				{
//					g_link_vector[li].TollMAP[dt] = 0;
//				}
//
//			}  //end of pt
//		}
//	}
//}


void g_ReadDemandFileBasedOnDemandFileList(Assignment& assignment)
{

	for (int i = 0; i < g_node_vector.size(); i++)
	{

		if (g_node_vector[i].zone_id >= 1 && assignment.g_zoneid_to_zone_seq_no_mapping.find(g_node_vector[i].zone_id) == assignment.g_zoneid_to_zone_seq_no_mapping.end())  // create a new zone  // we assume each zone only has one node
		{ // we need to make sure we only create a zone in the memory if only there is positive demand flow from the (new) OD table
			COZone ozone;
			ozone.node_seq_no = g_node_vector[i].node_seq_no;
			ozone.zone_id = g_node_vector[i].zone_id;
			ozone.zone_seq_no = g_zone_vector.size();
			assignment.g_zoneid_to_zone_seq_no_mapping[ozone.zone_id] = assignment.g_zoneid_to_zone_seq_no_mapping.size();  // create the zone id to zone seq no mapping

			g_zone_vector.push_back(ozone);  // add element into vector
											 //	cout << ozone.zone_id << ' ' << ozone.zone_seq_no << endl;
		}
	}

	cout << "number of zones = " << g_zone_vector.size() << endl;
	fprintf(g_pFileOutputLog, "number of zones =,%d\n", g_zone_vector.size());

	assignment.InitializeDemandMatrix(g_zone_vector.size(), assignment.g_DemandTypeVector.size());

	float total_demand_in_demand_file = 0;

	CCSVParser parser;
	cout << "Step 4: Reading file demand_file_list.csv..." << endl;

	if (parser.OpenCSVFile("demand_file_list.csv", true))
	{
		int i = 0;

		while (parser.ReadRecord())
		{

			int file_sequence_no = 1;
			string file_name;
			string format_type = "null";

			string demand_period, demand_type;

			int demand_format_flag = 0;

			if (parser.GetValueByFieldName("file_sequence_no", file_sequence_no) == false)
				break;

			if (file_sequence_no <= -1)  // skip negative sequence no 
				continue;

			parser.GetValueByFieldName("file_name", file_name);

			parser.GetValueByFieldName("demand_period", demand_period);



			parser.GetValueByFieldName("format_type", format_type);
			if (format_type.find("null") != string::npos)  // skip negative sequence no 
			{
				cout << "Please provide format_type in file demand_file_list.csv" << endl;
				g_ProgramStop();
			}



			double total_ratio = 0;

			parser.GetValueByFieldName("demand_type", demand_type);



			int demand_type_no = 0;
			int demand_period_no = 0;

			if (assignment.demand_period_to_seqno_mapping.find(demand_period) != assignment.demand_period_to_seqno_mapping.end())
			{
				demand_period_no = assignment.demand_period_to_seqno_mapping[demand_period];

			}

			if (assignment.demand_type_2_seqno_mapping.find(demand_period) != assignment.demand_type_2_seqno_mapping.end())
			{
				demand_type_no = assignment.demand_type_2_seqno_mapping[demand_period];

			}

			if (demand_period_no > _MAX_TIMEPERIODS)
			{
				cout << "demand_period_no should be less than settings in demand_period.csv. Please change the parameter settings in the source code." << endl;
				g_ProgramStop();
			}

			if (format_type.find("column") != string::npos)  // or muliti-column
			{


				bool bFileReady = false;
				int i;

				FILE* st;
				// read the file formaly after the test. 

				fopen_s(&st, file_name.c_str(), "r");
				if (st != NULL)
				{

					bFileReady = true;
					int line_no = 0;

					while (true)
					{
						int origin_zone = g_read_integer(st, true);

						if (origin_zone <= 0)
						{

							if (line_no == 1 && !feof(st))  // read only one line, but has not reached the end of the line
							{
								cout << endl << "Error: Only one line has been read from file. Are there multiple columns of demand type in file " << file_name << " per line?" << endl;
								g_ProgramStop();

							}
							break;
						}

						if (assignment.g_zoneid_to_zone_seq_no_mapping.find(origin_zone) == assignment.g_zoneid_to_zone_seq_no_mapping.end())
						{
							cout << endl << "Warning: origin zone " << origin_zone << "  has not been defined in input_node.csv" << endl;

							continue; // origin zone  has not been defined, skipped. 
						}

						int destination_zone = g_read_integer(st, true);

						if (assignment.g_zoneid_to_zone_seq_no_mapping.find(destination_zone) == assignment.g_zoneid_to_zone_seq_no_mapping.end())
						{
							cout << endl << "Warning: destination zone " << destination_zone << "  has not been defined in input_node.csv" << endl;

							continue; // destination zone  has not been defined, skipped. 
						}

						int from_zone_seq_no = 0;
						int to_zone_seq_no = 0;
						from_zone_seq_no = assignment.g_zoneid_to_zone_seq_no_mapping[origin_zone];
						to_zone_seq_no = assignment.g_zoneid_to_zone_seq_no_mapping[destination_zone];


						float demand_value = g_read_float(st);

						if (demand_value < -99) // encounter return 
						{
							break;
						}

						assignment.total_demand[demand_type_no][demand_period_no] += demand_value;
						assignment.g_demand_array[from_zone_seq_no][to_zone_seq_no][demand_type_no][demand_period_no].volume += demand_value;
						assignment.total_demand_volume += demand_value;
						assignment.g_origin_demand_array[from_zone_seq_no][demand_type_no][demand_period_no] += demand_value;

						// we generate vehicles here for each OD data line
						if (line_no <= 5)  // read only one line, but has not reached the end of the line
							cout << "origin:" << origin_zone << ", destination: " << destination_zone << ", value = " << demand_value << endl;

						if (line_no % 100000 == 0)
						{
							cout << "Reading file no." << file_sequence_no << ": " << file_name << " at " << line_no / 1000 << "K lines..." << endl;
						}


						line_no++;
					}  // scan lines


					fclose(st);

					cout << "total_demand_volume is " << assignment.total_demand_volume << endl;
				}
				else  //open file
				{
					cout << "Error: File " << file_name << " cannot be opened.\n It might be currently used and locked by EXCEL." << endl;
					g_ProgramStop();

				}
			}

			else if (format_type.compare("agent_csv") == 0)
			{
				vector<int> LineIntegerVector;

				CCSVParser parser;

				if (parser.OpenCSVFile(file_name, false))
				{
					int total_demand_in_demand_file = 0;

					int origin_zone, destination_zone;

					while (parser.ReadRecord())
					{
						total_demand_in_demand_file++;

						if (total_demand_in_demand_file % 1000 == 0)
							cout << "demand_volume is " << total_demand_in_demand_file << endl;

						parser.GetValueByFieldName("from_zone_id", origin_zone);
						parser.GetValueByFieldName("to_zone_id", destination_zone);

						if (assignment.g_zoneid_to_zone_seq_no_mapping.find(origin_zone) == assignment.g_zoneid_to_zone_seq_no_mapping.end())
						{
							continue; // origin zone  has not been defined, skipped. 
						}

						if (assignment.g_zoneid_to_zone_seq_no_mapping.find(destination_zone) == assignment.g_zoneid_to_zone_seq_no_mapping.end())
						{
							continue; // origin zone  has not been defined, skipped. 
						}

						float value = 1.0f;

						parser.GetValueByFieldName("volume", value);

						int from_zone_seq_no = 0;
						int to_zone_seq_no = 0;
						from_zone_seq_no = assignment.g_zoneid_to_zone_seq_no_mapping[origin_zone];
						to_zone_seq_no = assignment.g_zoneid_to_zone_seq_no_mapping[destination_zone];

						float number_of_vehicles = value;  // read the value

						int type = 1;  // first demand type definition
						int k_path_index = -1;
						parser.GetValueByFieldName("demand_type", type);
						parser.GetValueByFieldName("k_path_index", k_path_index);

						///// Jun: read agents' socio-demographic and other information
						int agent_id;
						int income;
						int gender;
						int flexibility;
						int vehicle;
						int purpose;
						float preferred_arrival_time;
						float actual_arrival_time;
						float travel_time_in_min;
						float free_flow_travel_time;

						parser.GetValueByFieldName("agent_id", agent_id);
						parser.GetValueByFieldName("income", income);
						parser.GetValueByFieldName("gender", gender);
						parser.GetValueByFieldName("flexibility", flexibility);
						parser.GetValueByFieldName("vehicle", vehicle);
						parser.GetValueByFieldName("purpose", purpose);
						parser.GetValueByFieldName("preferred_arrival_time", preferred_arrival_time);
						parser.GetValueByFieldName("travel_time_in_min", travel_time_in_min);
						parser.GetValueByFieldName("free_flow_travel_time", free_flow_travel_time);
						parser.GetValueByFieldName("arrival_time_in_min", actual_arrival_time);



						int demand_period_no = 0;
						parser.GetValueByFieldName("departure_time_period", demand_period_no);



						CAGBMAgent agent;

						agent.from_zone_seq_no = from_zone_seq_no;
						agent.to_zone_seq_no = to_zone_seq_no;
						agent.type = type;
						agent.time_period = demand_period_no;
						agent.k_path = k_path_index;
						agent.volume = number_of_vehicles;

						agent.agent_id = agent_id;
						agent.income = income;
						agent.flexibility = flexibility;
						agent.gender = gender;
						agent.vehicle = vehicle;
						agent.purpose = purpose;
						agent.preferred_arrival_time = preferred_arrival_time;
						agent.travel_time_in_min = travel_time_in_min;
						agent.arrival_time_in_min = actual_arrival_time;
						agent.free_flow_travel_time = free_flow_travel_time;

						g_agbmagent_vector.push_back(agent);

						assignment.total_demand[demand_type_no][demand_period_no] += number_of_vehicles;
						assignment.g_demand_array[from_zone_seq_no][to_zone_seq_no][demand_type_no][demand_period_no].volume += number_of_vehicles;
						assignment.total_demand_volume += number_of_vehicles;
						assignment.g_origin_demand_array[from_zone_seq_no][type][demand_period_no] += number_of_vehicles;
					}


				}
				else  //open file
				{
					cout << "Error: File " << file_name << " cannot be opened.\n It might be currently used and locked by EXCEL." << endl;
					g_ProgramStop();

				}

			}

			else
			{
				cout << "Error: format_type = " << format_type << " is not supported. Currently DTALite supports multi_column, matrix, full_matrix, dynasmart, agent_csv, agent_bin, trip_csv,transims_trip_file." << endl;
				g_ProgramStop();
			}
		}

	}

}



void g_ReadInputData(Assignment& assignment)
{

	//step 0:read demand period file
	CCSVParser parser_demand_period;
	cout << "Step 1: Reading file demand_period.csv..." << endl;
	//g_LogFile << "Step 7.1: Reading file input_demand_type.csv..." << g_GetUsedMemoryDataInMB() << endl;
	if (!parser_demand_period.OpenCSVFile("demand_period.csv", true))
	{
		cout << "demand_period.csv cannot be opened. " << endl;
		g_ProgramStop();

	}

	if (parser_demand_period.inFile.is_open() || parser_demand_period.OpenCSVFile("demand_period.csv", true))
	{

		while (parser_demand_period.ReadRecord())
		{

			CDemand_Period demand_period;


			if (parser_demand_period.GetValueByFieldName("demand_period_id", demand_period.demand_period_id) == false)
				break;

			if (parser_demand_period.GetValueByFieldName("demand_period", demand_period.demand_period) == false)
				break;


			vector<float> global_minute_vector;

			if (parser_demand_period.GetValueByFieldName("time_period", demand_period.time_period) == false)
				break;

			vector<string> input_string;
			input_string.push_back(demand_period.time_period);
			//input_string includes the start and end time of a time period with hhmm format
			global_minute_vector = g_time_parser(input_string); //global_minute_vector incldue the starting and ending time
			if (global_minute_vector.size() == 2)
			{

				demand_period.starting_time_slot_no = global_minute_vector[0] / MIN_PER_TIMESLOT;
				demand_period.ending_time_slot_no = global_minute_vector[1] / MIN_PER_TIMESLOT;

				//cout << global_minute_vector[0] << endl;
				//cout << global_minute_vector[1] << endl;
			}

			assignment.demand_period_to_seqno_mapping[demand_period.demand_period] = assignment.g_DemandPeriodVector.size();

			assignment.g_DemandPeriodVector.push_back(demand_period);


		}
		parser_demand_period.CloseCSVFile();
	}
	else
	{
		cout << "Error: File input_demand_type.csv cannot be opened.\n It might be currently used and locked by EXCEL." << endl;
		g_ProgramStop();
	}


	assignment.g_number_of_demand_periods = assignment.g_DemandPeriodVector.size();
	//step 1:read demand type file
	CCSVParser parser_demand_type;
	cout << "Step 2: Reading file demand_type.csv..." << endl;
	//g_LogFile << "Step 7.1: Reading file input_demand_type.csv..." << g_GetUsedMemoryDataInMB() << endl;
	if (!parser_demand_type.OpenCSVFile("demand_type.csv", true))
	{
		cout << "demand_type.csv cannot be opened. " << endl;
		g_ProgramStop();

	}

	if (parser_demand_type.inFile.is_open() || parser_demand_type.OpenCSVFile("demand_type.csv", true))
	{
		assignment.g_DemandTypeVector.clear();
		while (parser_demand_type.ReadRecord())
		{

			CDemand_Type demand_type;

			if (parser_demand_type.GetValueByFieldName("demand_type_id", demand_type.demand_type_id) == false)
				break;

			if (parser_demand_type.GetValueByFieldName("demand_type", demand_type.demand_type) == false)
				break;

			parser_demand_type.GetValueByFieldName("VOT", demand_type.value_of_time);
			parser_demand_type.GetValueByFieldName("cap_ratio", demand_type.capacity_equivalent_ratio);


			assignment.demand_type_2_seqno_mapping[demand_type.demand_type] = assignment.g_DemandTypeVector.size();



			assignment.g_DemandTypeVector.push_back(demand_type);
			assignment.g_number_of_demand_types = assignment.g_DemandTypeVector.size();

		}
		parser_demand_type.CloseCSVFile();
	}
	else
	{
		cout << "Error: File input_demand_type.csv cannot be opened.\n It might be currently used and locked by EXCEL." << endl;
		g_ProgramStop();
	}


	if (assignment.g_DemandTypeVector.size() >= _MAX_DEMANDTYPES)
	{
		cout << "Error: demand_type = " << assignment.g_DemandTypeVector.size() << " in file demand_type.csv is too large. " << "_MAX_DEMANDTYPES = " << _MAX_DEMANDTYPES << "Please contact program developers!";

		g_ProgramStop();
	}


	assignment.g_number_of_nodes = 0;
	assignment.g_number_of_links = 0;  // initialize  the counter to 0

	if (assignment.g_number_of_K_paths > _MAX_K_PATH)
	{
		cout << "g_number_of_K_paths >= _MAX_K_PATH" << endl;
		g_ProgramStop();
	}

	int internal_node_seq_no = 0;
	double x, y;
	// step 3: read node file 
	CCSVParser parser;
	if (parser.OpenCSVFile("node.csv", true))
	{

		while (parser.ReadRecord())  // if this line contains [] mark, then we will also read field headers.
		{

			string node_id;

			if (parser.GetValueByFieldName("node_id", node_id) == false)
				continue;

			if (assignment.g_internal_node_to_seq_no_map.find(node_id) != assignment.g_internal_node_to_seq_no_map.end())
			{
				continue; //has been defined
			}
			assignment.g_internal_node_to_seq_no_map[node_id] = internal_node_seq_no;


			CNode node;  // create a node object 

			node.node_id = node_id;
			node.node_seq_no = internal_node_seq_no;
			parser.GetValueByFieldName("zone_id", node.zone_id);



			/*node.x = x;
			node.y = y;*/
			internal_node_seq_no++;

			g_node_vector.push_back(node);  // push it to the global node vector

			assignment.g_number_of_nodes++;
			if (assignment.g_number_of_nodes % 1000 == 0)
				cout << "reading " << assignment.g_number_of_nodes << " nodes.. " << endl;
		}

		cout << "number of nodes = " << assignment.g_number_of_nodes << endl;

		fprintf(g_pFileOutputLog, "number of nodes =,%d\n", assignment.g_number_of_nodes);


		parser.CloseCSVFile();
	}


	// step 4: read link file 

	CCSVParser parser_link;

	if (parser_link.OpenCSVFile("road_link.csv", true))
	{
		while (parser_link.ReadRecord())  // if this line contains [] mark, then we will also read field headers.
		{
			string from_node_id;
			string to_node_id;
			if (parser_link.GetValueByFieldName("from_node_id", from_node_id) == false)
				continue;
			if (parser_link.GetValueByFieldName("to_node_id", to_node_id) == false)
				continue;

			string linkID = "0";
			parser_link.GetValueByFieldName("road_link_id", linkID);


			// add the to node id into the outbound (adjacent) node list

			int internal_from_node_seq_no = assignment.g_internal_node_to_seq_no_map[from_node_id];  // map external node number to internal node seq no. 
			int internal_to_node_seq_no = assignment.g_internal_node_to_seq_no_map[to_node_id];

			CLink link;  // create a link object 

			link.from_node_seq_no = internal_from_node_seq_no;
			link.to_node_seq_no = internal_to_node_seq_no;
			link.link_seq_no = assignment.g_number_of_links;
			link.to_node_seq_no = internal_to_node_seq_no;
			link.link_id = linkID;

			parser_link.GetValueByFieldName("facility_type", link.type);

			string demand_type_code;
			if (!parser_link.GetValueByFieldName("demand_type_code", demand_type_code))
				demand_type_code = "";
			link.demand_type_code = demand_type_code;

			float length = 1.0; // km or mile
			float free_speed = 1.0;
			float k_jam = 200;
			parser_link.GetValueByFieldName("length", length);
			parser_link.GetValueByFieldName("free_speed", free_speed);

			int number_of_lanes = 1;
			parser_link.GetValueByFieldName("lanes", number_of_lanes);

			char VDF_field_name[20];

			for (int tau = 0; tau < assignment.g_number_of_demand_periods; tau++)
			{
				int demand_period_id = assignment.g_DemandPeriodVector[tau].demand_period_id;
				sprintf_s (VDF_field_name, "VDF_fftt%d", demand_period_id);
				parser_link.GetValueByFieldName(VDF_field_name, link.VDF_period[tau].FFTT);

				sprintf_s (VDF_field_name, "VDF_cap%d", demand_period_id);
				parser_link.GetValueByFieldName(VDF_field_name, link.VDF_period[tau].capacity);

				sprintf_s (VDF_field_name, "VDF_alpha%d", demand_period_id);
				parser_link.GetValueByFieldName(VDF_field_name, link.VDF_period[tau].alpha);

				sprintf_s (VDF_field_name, "VDF_beta%d", demand_period_id);
				parser_link.GetValueByFieldName(VDF_field_name, link.VDF_period[tau].beta);


				sprintf_s (VDF_field_name, "VDF_mu%d", demand_period_id);
				parser_link.GetValueByFieldName(VDF_field_name, link.VDF_period[tau].mu);

				sprintf_s (VDF_field_name, "VDF_gamma%d", demand_period_id);
				parser_link.GetValueByFieldName(VDF_field_name, link.VDF_period[tau].gamma);

				link.VDF_period[tau].starting_time_slot_no = assignment.g_DemandPeriodVector[tau].starting_time_slot_no;
				link.VDF_period[tau].ending_time_slot_no = assignment.g_DemandPeriodVector[tau].ending_time_slot_no;
				link.VDF_period[tau].period = assignment.g_DemandPeriodVector[tau].time_period;


			}
			// for each period


			float default_cap = 1000;
			float default_BaseTT = 1;

			// setup default value
			for (int tau = 0; tau < assignment.g_number_of_demand_periods; tau++)
			{
				link.TDBaseTT[tau] = default_BaseTT;
				link.TDBaseCap[tau] = default_cap;
			}

			//link.m_OutflowNumLanes = number_of_lanes;//visum lane_cap is actually link_cap

			link.link_spatial_capacity = k_jam * number_of_lanes*length;

			link.free_flow_travel_time_in_min = default_BaseTT;

			link.length = length;
			for (int tau = 0; tau < assignment.g_number_of_demand_periods; tau++)
			{
				link.travel_time_per_period[tau] = length / free_speed * 60;
			}
			// min // calculate link cost based length and speed limit // later we should also read link_capacity, calculate delay 


			g_node_vector[internal_from_node_seq_no].m_outgoing_link_vector.push_back(link);  // add this link to the corresponding node as part of outgoing node/link


			g_link_vector.push_back(link);

			assignment.g_number_of_links++;

			if (assignment.g_number_of_links % 1000 == 0)
				cout << "reading " << assignment.g_number_of_links << " links.. " << endl;
		}

		parser_link.CloseCSVFile();
	}
	// we now know the number of links
	cout << "number of links = " << assignment.g_number_of_links << endl;

	fprintf(g_pFileOutputLog, "number of links =,%d\n", assignment.g_number_of_links);

	parser_link.CloseCSVFile();


}


float total_tree_cost[_MAX_K_PATH];
float total_tree_distance[_MAX_K_PATH];
float total_pi_cost[_MAX_K_PATH];
float total_experienced_cost[_MAX_K_PATH];

float _gap_[_MAX_K_PATH];
float _gap_relative_[_MAX_K_PATH];



void g_output_simulation_result(Assignment& assignment)
{
	//int b_debug_detail_flag = 0;
	//FILE* g_pFileLinkMOE = NULL;
	//g_pFileLinkMOE = fopen_s("link_performance.csv", "w");
	//if (g_pFileLinkMOE == NULL)
	//{
	//	cout << "File link_performance.csv cannot be opened." << endl;
	//	g_ProgramStop();
	//}
	//else
	//{
	//	// Option 1: BPR function
	//	fprintf(g_pFileLinkMOE, "link_id,time_period,volume,travel_time,density,voc_ratio,notes\n");

	//	for (int l = 0; l < g_link_vector.size(); l++) //Initialization for all nodes
	//	{
	//		for (int tau = 0; tau < assignment.g_number_of_demand_periods; tau++)
	//		{
	//			CLink* pLink = &(g_link_vector[l]);

	//			fprintf(g_pFileLinkMOE, "%s,%s,%.3f,%.3f,%.1f,%.1f,%s\n",
	//				pLink->link_id.c_str(),
	//				assignment.g_DemandPeriodVector[tau].demand_period.c_str(),
	//				pLink->flow_volume_per_period[tau],
	//				pLink->travel_time_per_period[tau],
	//				-1,
	//				pLink->get_VOC_ratio(tau),
	//				"BPR model");

	//		}
	//	}
	//}
	//fclose(g_pFileLinkMOE);

	int b_debug_detail_flag = 0;
	FILE* g_pFileLinkMOE = NULL;
	fopen_s(&g_pFileLinkMOE,"link_performance.csv", "w");
	if (g_pFileLinkMOE == NULL)
	{
		cout << "File link_performance.csv cannot be opened." << endl;
		g_ProgramStop();
	}
	else
	{
		// Option 2: BPR_X function
		fprintf(g_pFileLinkMOE, "road_link_id,time_period,volume,travel_time,notes\n");
		for (int l = 0; l < g_link_vector.size(); l++) //Initialization for all nodes
		{
			for (int tau = 0; tau < assignment.g_number_of_demand_periods; tau++)
			{
				int ii = FLAG_PEORIOD_TIMESLOT;

				if (ii == 0)
				{
					fprintf(g_pFileLinkMOE, "%s,%s,%.3f,%.3f,%s\n",
						g_link_vector[l].link_id.c_str(),
						g_link_vector[l].VDF_period[tau].period.c_str(),
						g_link_vector[l].flow_volume_per_period[tau],
						g_link_vector[l].VDF_period[tau].avg_travel_time,
						"period-based");
				}
				else
				{
					int starting_time = g_link_vector[l].VDF_period[tau].starting_time_slot_no;
					int ending_time = g_link_vector[l].VDF_period[tau].ending_time_slot_no;

					for (int t = 0; t <= ending_time - starting_time; t++)
					{
						fprintf(g_pFileLinkMOE, "%s,%s,%.3f,%.3f,%s\n",
							g_link_vector[l].link_id.c_str(),
							assignment.g_DemandPeriodVector[tau].demand_period.c_str(),
							g_link_vector[l].VDF_period[tau].arrival_rate[t],
//							g_link_vector[l].VDF_period[tau].waiting_time[t] + g_link_vector[l].VDF_period[tau].FFTT,
							g_link_vector[l].VDF_period[tau].get_waiting_time(t) + g_link_vector[l].VDF_period[tau].FFTT,

							
							"timeslot-dependent");
					}

				}

			}

		}
	}

	fclose(g_pFileLinkMOE);

	//added by zhuge,to output the ODMOE 
	FILE* g_pFileODMOE = NULL;
	fopen_s(&g_pFileODMOE,"agent.csv", "w");
	if (g_pFileODMOE == NULL)
	{
		cout << "File od_performance.csv cannot be opened." << endl;
		g_ProgramStop();
	}
	else
	{
		fprintf(g_pFileODMOE, "agent_id,from_zone_id,to_zone_id,from_node_id,to_node_id,demand_type,demand_period,volume,cost,travel_time,distance,node_sequence\n");

		int count = 1;
		for (int i = 0; i < g_zone_vector.size(); i++)
			for (int j = 0; j < g_zone_vector.size(); j++)
				for (int dt = 0; dt < assignment.g_DemandTypeVector.size(); dt++)
					for (int tau = 0; tau < assignment.g_DemandPeriodVector.size(); tau++)
					{

						if (assignment.g_demand_array[i][j][dt][tau].volume > 0)
						{
							fprintf(g_pFileODMOE, "%d,%d,%d,%s,%s,%s,%s,%.1f,%.1f,%.1f,%.4f,",
								count,
								g_zone_vector[i].zone_id,
								g_zone_vector[j].zone_id,
								g_node_vector[g_zone_vector[i].node_seq_no].node_id.c_str(),
								g_node_vector[g_zone_vector[j].node_seq_no].node_id.c_str(),
								assignment.g_DemandTypeVector[dt].demand_type.c_str(),
								assignment.g_DemandPeriodVector[tau].demand_period.c_str(),
								assignment.g_demand_array[i][j][dt][tau].volume,
								assignment.g_demand_array[i][j][dt][tau].cost,
								assignment.g_demand_array[i][j][dt][tau].time,
								assignment.g_demand_array[i][j][dt][tau].distance

							);


							for (int ni = assignment.g_demand_array[i][j][dt][tau].path_node_vector.size() - 1; ni >= 0; ni--)
							{
								fprintf(g_pFileODMOE, "%s;", g_node_vector[assignment.g_demand_array[i][j][dt][tau].path_node_vector[ni]].node_id.c_str());

							}

							count++;
							fprintf(g_pFileODMOE, "\n");
						}
					}
	}
	fclose(g_pFileODMOE);
}

//***
// major function 1:  allocate memory and initialize the data
// void AllocateMemory(int number_of_nodes)
//
//major function 2: // time-dependent label correcting algorithm with double queue implementation
//int optimal_label_correcting(int origin_node, int destination_node, int departure_time, int shortest_path_debugging_flag, FILE* pFileDebugLog, Assignment& assignment, int time_period_no = 0, int demand_type = 1, float VOT = 10)

//	//major function: update the cost for each node at each SP tree, using a stack from the origin structure 
//int tree_cost_updating(int origin_node, int departure_time, int shortest_path_debugging_flag, FILE* pFileDebugLog, Assignment& assignment, int time_period_no = 0, int demand_type = 1)

//***

// The one and only application object
using namespace std;

int g_number_of_CPU_threads()
{
	int number_of_threads = omp_get_max_threads();

	int max_number_of_threads = 4000;

	if (number_of_threads > max_number_of_threads)
		number_of_threads = max_number_of_threads;

	return number_of_threads;

}


void g_assign_computing_tasks_to_memory_blocks(Assignment& assignment)
{
	//fprintf(g_pFileDebugLog, "-------g_assign_computing_tasks_to_memory_blocks-------\n");
	// step 2: assign node to thread
	for (int k = 0; k < assignment.g_number_of_K_paths; k++)
	{
		for (int i = 0; i < g_node_vector.size(); i++)  //assign all nodes to the corresponding thread
		{
			if (g_node_vector[i].zone_id >= 1)
			{
				int zone_seq_no = assignment.g_zoneid_to_zone_seq_no_mapping[g_node_vector[i].zone_id];


				for (int dt = 0; dt < assignment.g_DemandTypeVector.size(); dt++)
				{

					for (int tau = 0; tau < assignment.g_DemandPeriodVector.size(); tau++)
					{
						if (assignment.g_origin_demand_array[zone_seq_no][dt][tau] > 0) // with feasible flow
						{
							//fprintf(g_pFileDebugLog, "%f\n",g_origin_demand_array[zone_seq_no][dt][tau]);

								//cout << assignment.g_origin_demand_array[zone_seq_no][dt][tau] << endl;

							NetworkForSP* p_NetworkForSP = new NetworkForSP();
							g_NetworkForSP_vector.push_back(p_NetworkForSP);
							p_NetworkForSP->m_iteration_k = k;
							p_NetworkForSP->m_origin_node = i;
							p_NetworkForSP->m_origin_zone_seq_no = zone_seq_no;

							p_NetworkForSP->m_demand_type_no = dt;
							p_NetworkForSP->m_demand_time_period_no = tau;
							p_NetworkForSP->AllocateMemory(assignment.g_number_of_nodes);

						}
					}
				}


			}
		}
	}

}



void g_reset_link_volume(int number_of_links)
{
	for (int l = 0; l < number_of_links; l++)
	{
		for (int tau = 0; tau < assignment.g_number_of_demand_periods; tau++)
		{
			g_link_vector[l].flow_volume_per_period[tau] = 0;
			g_link_vector[l].queue_length_perslot[tau] = 0;
		}
	}

}

//major function: update the cost for each node at each SP tree, using a stack from the origin structure 

int NetworkForSP::calculate_TD_link_flow(Assignment& assignment, int iteration_number_outterloop)
{

	int origin_node = m_origin_node;
	int departure_time = m_demand_time_period_no;
	int shortest_path_debugging_flag = 0;
	int demand_type = m_demand_type_no;


	if (m_iteration_k > iteration_number_outterloop)  // we only update available path tree;
		return 0;

	if (g_node_vector[origin_node].m_outgoing_link_vector.size() == 0)
	{
		return 0;
	}

	//	fprintf(g_pFileDebugLog, "------------START: origin:  %d  ; Departure time: %d  ; demand type:  %d  --------------\n", origin_node + 1, departure_time, demand_type);
	float k_path_prob = float(1) / float(iteration_number_outterloop + 1);  //XZ: use default value as MSA
	int num = 0;
	for (int i = 0;i < assignment.g_number_of_nodes;i++)
	{

		if (g_node_vector[i].zone_id >= 1)
		{
			//			fprintf(g_pFileDebugLog, "--------origin  %d ; destination node: %d ; (zone: %d) -------\n", origin_node + 1, i+1, g_node_vector[i].zone_id);
			//fprintf(g_pFileDebugLog, "--------iteration number outterloop  %d ;  -------\n", iteration_number_outterloop);
			int destination_zone_seq_no = assignment.g_zoneid_to_zone_seq_no_mapping[g_node_vector[i].zone_id];

			float car_ratio = 1;// assignment.g_DemandTypeVector[demand_type].capacity_equivalent_ratio;
			float volume = assignment.g_demand_array[m_origin_zone_seq_no][destination_zone_seq_no][demand_type][m_demand_time_period_no].volume * car_ratio*k_path_prob;

			if (volume > 0.000001)
			{
				vector<int> temp_path_node_vector; //node seq vector for each ODK, zhuge

				int current_node_seq_no = i;  // destination node
				int current_link_seq_no = -1;

				// backtrace the sp tree from the destination to the root (at origin) 
				while (current_node_seq_no >= 0 && current_node_seq_no < assignment.g_number_of_nodes)
				{

					temp_path_node_vector.push_back(current_node_seq_no);
					if (current_node_seq_no >= 0 && current_node_seq_no < assignment.g_number_of_nodes)  // this is valid node 
					{
						current_link_seq_no = m_link_predecessor[current_node_seq_no];

						// fetch m_link_predecessor to mark the link volume

						if (current_link_seq_no >= 0 && current_link_seq_no < assignment.g_number_of_links)
						{
							//static
							//link_volume[current_link_seq_no] += (destination_k_path_prob[d] * destination_OD_volume[d][assignment.g_DemandTypeVector[i]][tau]);

							float value_before = g_link_vector[current_link_seq_no].flow_volume_per_period[m_demand_time_period_no];

							g_link_vector[current_link_seq_no].flow_volume_per_period[m_demand_time_period_no] += volume;


							////								fprintf(g_pFileDebugLog, "\n o_k = %d; t_k=%d; origin  %d ; destination node: %d ; (zone: %d); link (%d,%d, t=%d)---volume = %f + %f = %f", 
							//									iteration_number_outterloop, m_iteration_k, origin_node + 1, i + 1, g_node_vector[i].zone_id,
							//									g_link_vector[current_link_seq_no].from_node_seq_no + 1, g_link_vector[current_link_seq_no].to_node_seq_no + 1, m_demand_time_period_no,
							//									value_before, volume, g_link_vector[current_link_seq_no].flow_volume_per_period[m_demand_time_period_no]);
						}
					}
					current_node_seq_no = m_node_predecessor[current_node_seq_no];  // update node seq no	
				}
				//fprintf(g_pFileDebugLog, "\n");

				// we obtain the cost, time, distance from the last tree-k 
				if (m_iteration_k == assignment.g_number_of_K_paths - 1)  // last tree-k iteration for the total number of iterations
				{
					assignment.g_demand_array[m_origin_zone_seq_no][destination_zone_seq_no][demand_type][m_demand_time_period_no].cost = m_node_label_cost[i];
					assignment.g_demand_array[m_origin_zone_seq_no][destination_zone_seq_no][demand_type][m_demand_time_period_no].time = m_label_time_array[i];
					assignment.g_demand_array[m_origin_zone_seq_no][destination_zone_seq_no][demand_type][m_demand_time_period_no].distance = m_label_distance_array[i];
					assignment.g_demand_array[m_origin_zone_seq_no][destination_zone_seq_no][demand_type][m_demand_time_period_no].path_node_vector = temp_path_node_vector;
				}

			}
		}

	}

	//fprintf(g_pFileDebugLog, "-------------END---------------- \n");
}



void  CLink::CalculateTD_VDFunction()
{


	for (int tau = 0; tau < assignment.g_number_of_demand_periods; tau++)
	{
		float starting_time_slot_no = assignment.g_DemandPeriodVector[tau].starting_time_slot_no;
		float ending_time_slot_no = assignment.g_DemandPeriodVector[tau].ending_time_slot_no;
		travel_time_per_period[tau] = VDF_period[tau].PeroformBPR(flow_volume_per_period[tau]);
		//travel_time_per_period[tau] = VDF_period[tau].PeroformBPR_X(flow_volume_per_period[tau]);

	}
}

double network_assignment(int iteration_number, int b)
{

	fopen_s(&g_pFileOutputLog, "output_solution.csv", "w");
	if (g_pFileOutputLog == NULL)
	{
		cout << "File output_solution.csv cannot be opened." << endl;
		g_ProgramStop();
	}

	assignment.g_number_of_K_paths = iteration_number;


	// step 1: read input data of network / demand tables / Toll
	g_ReadInputData(assignment);
	g_ReadDemandFileBasedOnDemandFileList(assignment);

	// definte timestamps
	clock_t start_t, end_t, total_t, start_t_1, end_t_1, total_t_1;
	int i;

	//step 2: allocate memory and assign computing tasks
	g_assign_computing_tasks_to_memory_blocks(assignment);

	//step 3: find shortest path and update path cost of tree using TD travel time
	for (int iteration_number = 0; iteration_number < assignment.g_number_of_K_paths; iteration_number++)
	{
		//TRACE("Loop 1: assignment iteration %d", iteration_number);
		//step 3: compute K SP tree
		start_t = clock();

		g_reset_link_volume(g_link_vector.size());  // we can have a recursive formulat to reupdate the current link volume by a factor of k/(k+1), and use the newly generated path flow to add the additional 1/(k+1)
		//g_setup_link_cost_in_each_memory_block(iteration_number, assignment);
#pragma omp parallel for  // step 3: C++ open mp automatically create n threads., each thread has its own computing thread on a cpu core 
		for (int ProcessID = 0; ProcessID < g_NetworkForSP_vector.size(); ProcessID++) //changed by zhuge
		{
			g_NetworkForSP_vector[ProcessID]->optimal_label_correcting(assignment, iteration_number);
		}


		for (int ProcessID = 0; ProcessID < g_NetworkForSP_vector.size(); ProcessID++)
		{
			g_NetworkForSP_vector[ProcessID]->calculate_TD_link_flow(assignment, iteration_number);
		}


		end_t = clock();
		total_t = (end_t - start_t);

		//					cout << "CPU Running Time for SP = " << total_t << " milliseconds" << " " << "k=" << iteration_number << endl;
		//					fprintf(g_pFileDebugLog, "CPU Running Time for SP  = %ld milliseconds\n", total_t);

		////step 3.2: calculate TD link travel time using TD inflow flow and capacity  
		//					start_t_1 = clock();

#pragma omp parallel for  // step collect all partial link volume to compute link volume across all zones
		for (int l = 0; l < g_link_vector.size(); l++)
		{
			//g_link_vector[l].tally_flow_volume_across_all_processors();
			g_link_vector[l].CalculateTD_VDFunction();

		}

	}
	cout << "Output for old demand results Done!" << endl;

	cout << "CPU Running Time for assignment= " << total_t << " milliseconds" << endl;
	//step 4: output simulation results of the new demand 
	g_output_simulation_result(assignment);

	cout << "CPU Running Time for Reassignment= " << total_t << " milliseconds" << endl;
	cout << "End of Optimization " << endl;
	cout << "free memory.." << endl;
	cout << "done." << endl;

	g_node_vector.clear();

	for (int l = 0; l < g_link_vector.size(); l++)
	{
		g_link_vector[l].free_memory();
	}
	g_link_vector.clear();
	getchar();
	return 1.0;

}
