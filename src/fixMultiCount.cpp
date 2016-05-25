#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string>
#include <vector>
#include <unordered_map>
#include <Rcpp.h>
#include <math.h>

using namespace std;
using namespace Rcpp;
typedef unordered_map<string, int> countDict;
typedef unordered_map<string, vector<string> > multiCountDict;
typedef vector<string> stringlist;
typedef vector<int> intlist;

intlist getCounts(stringlist gene_ids, countDict geneCountDict)
{
	// given a list of gene ids return the counts of these genes
	// base on the given count dict
	string geneID;
	intlist countList(gene_ids.size());
	for (int i = 0; i < gene_ids.size(); i++)
	{
		geneID = gene_ids[i];
		if (geneCountDict.find(geneID) == geneCountDict.end())
		{
			countList[i] = 0;
		}
		else
		{
			countList[i] = geneCountDict[geneID];
		}
	}
	return countList;
}

DataFrame constructResult(countDict geneCountDict)
{
	// convert count hash table to dataframe
	int geneCounts = geneCountDict.size();
	int i = 0;
	intlist counts(geneCounts);
	stringlist gene_ids(geneCounts);
	for( const auto& record : geneCountDict)
	{
		gene_ids[i] = record.first;
		counts[i] = record.second;
		i ++;
	}
	DataFrame newDF = DataFrame::create( Named("id")=gene_ids,
										Named("count") = counts);
	return newDF;
}


countDict updateCounts(countDict geneCountDict, multiCountDict multiCountsID)
// getting the list of genes that a single fragment mapped to
// and select the highest-count gene as assignment
{
	stringlist gene_ids;
	intlist::iterator maxIter;
	string addingGeneID;
	int position;
	int i = 0;
	for( const auto& record : multiCountsID )
	{
		gene_ids = record.second;
		intlist countList(gene_ids.size());
 		countList = getCounts(gene_ids, geneCountDict);
		maxIter = max_element(countList.begin(), countList.end());
		position = distance(countList.begin(), maxIter);
		addingGeneID = gene_ids[position];
	//	cout << addingGeneID << '\n';
		if (geneCountDict.find(addingGeneID) == geneCountDict.end())
		{
			geneCountDict[addingGeneID] = 1;
		}
		else
		{
			geneCountDict[addingGeneID] ++;
		}
		i ++;
		if (remainder(i,10000) == 0)
		{
			cout << "Assigned " << i << " multi-mapped fragments" << '\n';
		}
	}
	return geneCountDict;
}

int columnExist(stringlist neededName, stringlist existingNames, string dfname)
{
  string colname;
  string errorMessage;
  for (int i=0;i<neededName.size();i++)
  {
    colname = neededName[i];
    if (find(existingNames.begin(),existingNames.end(),colname) == existingNames.end())
    {
      errorMessage = "column ["+ colname +"] not in data frame:" + dfname + "\n";
      Rcpp::stop(errorMessage);
    }
  }
  return 0;
}

//' @useDynLib correctMultiCount
//' @importFrom Rcpp evalCpp
//' @exportPattern "^[[:alpha:]]+"
//' @title correctCounts
//' @description Assigning multiply-mapped reads to most-abundantly-mapped locus
//' @param basecount: should be a dataframe for storing gene count from uniquely-mapped reads **Have to contain columns named as 'id' for gene_id and 'count' for counts
//' @param   multiCount should be a dataframe for storing multiply-mapped reads and their mapped locus. Columns named as 'fragment_id' for storing read id  and 'gene_id' for storing gene_id
//' @export
// [[Rcpp::export]]
DataFrame correctCounts(DataFrame baseCount, DataFrame multiCount)
{
  //Check columns
  stringlist baseCol = baseCount.names(), multiCol = multiCount.names();
  stringlist baseColnames(2), multiColnames(2);
  baseColnames[0] = "id";
  baseColnames[1] = "count";
  multiColnames[0] = "fragment_id";
  multiColnames[1] = "gene_id";
  columnExist(baseColnames, baseCol, "baseCount");
  columnExist(multiColnames, multiCol, "multiCount");
	// read base count file and generate count dict
	countDict geneCountDict;
	IntegerVector counts = baseCount["count"];
	CharacterVector geneIDs = baseCount["id"];
	string id;
	int count;

	// creating hash table for gene counts
	for (int i = 0; i < counts.size();i++)
	{
	  id = geneIDs[i];
	  count = counts[i];
	  geneCountDict[id] = count;
	}
	cout << "Read gene count table with "<< geneCountDict.size() << " genes" << '\n';

	//----------------------------------------------------------------------------//
	// read multiple counts
	CharacterVector fragment_ids = multiCount["fragment_id"];
	CharacterVector multi_gene_ids = multiCount["gene_id"];
	string geneID;
	string fragment_ID;
	multiCountDict multiCountsID;
	for (int i = 0; i < fragment_ids.size(); i++)
	{
		geneID = multi_gene_ids[i];
		fragment_ID = fragment_ids[i];
		multiCountsID[fragment_ID].push_back(geneID);
	}
	cout << "Parsed " << multiCountsID.size() << " multi-mapped fragments" << '\n';

	// -----------------------------------------------------------------------------//
	// update count dict by getting the maximum-counted-gene
	geneCountDict = updateCounts(geneCountDict, multiCountsID);

	// -----------------------------------------------------------------------------//
	//
	// Create new count table
	DataFrame newDF = constructResult(geneCountDict);
	return newDF;
}

