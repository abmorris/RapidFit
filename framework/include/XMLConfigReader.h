/*!
 * @class XMLConfigReader
 *
 * @brief Opens an xml config file and uses it to create RapidFit data objects
 *
 * This can be thought of as the initialization engine of RapidFit
 *
 * It is capable of calling the correct constructors and assembling the correct objects for passing back into the main code
 *
 *
 * @ warning XMLConfigurationReader DOES NOT UNDER ANY CIRCUMSTANCES LOOK AFTER THE OBJECTS IT PASSES A POINTER FOR
 *                                  clean up after yourself!
 *
 *
 * This does rely on XMLTag but in principle it can be exchanged for a proper XML parser
 *
 * Benefits?	Unknown, none of us are XML guru's
 *
 * Cons?        Can't comment out code as quickly/flexibly
 *              Would be difficult to engineer --OverrideXML into a new parser it barely works here ;)
 *              Would have to be installed on all platforms and the authors have a distrust of ROOT XML
 *
 * @author Benjamin M Wynne bwynne@cern.ch
 */

#pragma once
#ifndef XML_CONFIG_READER_H
#define XML_CONFIG_READER_H

//	RapidFit Headers
#include "XMLTag.h"
#include "ParameterSet.h"
#include "PDFWithData.h"
#include "FitFunctionConfiguration.h"
#include "MinimiserConfiguration.h"
#include "OutputConfiguration.h"
#include "DataSetConfiguration.h"
#include "ComponentPlotter.h"
#include "ConstraintFunction.h"
#include "ScanParam.h"
#include "PrecalculatorConfig.h"
#include "DebugClass.h"
#include "I_XMLConfigReader.h"
//	System Headers
#include <vector>
#include <string>

using std::string;
using std::pair;
using std::vector;

class XMLConfigReader : public I_XMLConfigReader
{
	public:

		static void PrintTag( const XMLTag* thisTag );

		/*!
		 * @brief Constructor with Input FileName and optional override flags
		 *
		 * @param FileName   This is the name of the input file to parse.
		 *                   This will exit cleanly if none is found
		 *
		 * @param Override   This is a very expert feature of RapidFit
		 *                   This allows us to construct sequential tests modify some part of the XML file each run
		 *
		 */
		XMLConfigReader( string FileName, vector<pair<string,string> >* Override = NULL, string GlobalPrepend=string() );

		/*!
		 * Default Destructor
		 */
		virtual ~XMLConfigReader();

		/*!
		 * @brief Returns the result of a boolean decision to determine if the XML has passed some basic sanity checks for having a complete configuration
		 *
		 * @return True if the file can be used by RapidFit, False if it cannot
		 */
		bool IsValid() const;

		/*!
		 * @brief Return the Input XML from the file
		 *
		 * @return This has removed all whitespaces from the XML so output should be treated as a template
		 */
		vector<string> GetXML() const;


		/*!
		 * @brief Returns the ParameterSet for the whole Fit.
		 *
		 * When an XML File contains multiple ParameterSets it will attempt to merge them into one large ParameterSet
		 * with one of each unique Parameter
		 *
		 * When Parameters are Provided at the Command Line they will Override the existing Parameters in the XML or will Add to it
		 *
		 * @param CommandLineParam  This parameter contains strings from Runtime to override or add a Physics Parameter to the ParameterSet
		 *
		 * @return Returns a Pointer to the ParameterSet object and gives up any control over it.
		 */
		ParameterSet* GetFitParameters( vector<string> CommandLineParam = vector<string>() );

		/*!
		 * @brief Returns the MinimiserConfiguration class as defined in the XML
		 *
		 * @return Returns a Pointer to the MinimiserConfiguration class which is capable of Creating the Minimiser required for the Fit
		 */
		MinimiserConfiguration * GetMinimiserConfiguration();

		/*!
		 * @brief Returns the FitFunctionConfiguration class as defined in the XML
		 *
		 * @return Returns a Pointer to the FitFunctionConfiguration class which is capable of Creating the FitFunction required for the Fit
		 */
		FitFunctionConfiguration * GetFitFunctionConfiguration();

		/*!
		 * @brief Returns the OutputConfiguration as defined in the XML
		 *
		 * @return Returns a Pointer to the OutputConfiguration which produces more detailed output around the fit minima
		 */
		OutputConfiguration * GetOutputConfiguration();

		/*!
		 * @brief Returns a vector of PDFWithData objects, one for each ToFit containing a PDF in the XML
		 *
		 * @param StartingValues   This contains the starting values passed from the command line to control where in the file the fits read in data
		 *
		 * @return Returns pointers to the different PDFWithData objects present in the XML
		 */
		vector< PDFWithData* > GetPDFsAndData( vector<int> StartingValues = vector<int>() );

		/*!
		 * @brief Returns the list of External Constrains as defined in the XML
		 *
		 * @returns a vector of separate ConstraintFunctions as defined in the XML
		 */
		vector< ConstraintFunction* > GetConstraints();

		/*!
		 * @brief Returns a vector of the PhaseSpaceBoundaries associated with the PDFs in the XML
		 */
		vector< PhaseSpaceBoundary* > GetPhaseSpaceBoundaries();

		/*!
		 * @brief Construct the PrecalculatorConfig Object in order to Weight a DataSet using an algorithm
		 */
		PrecalculatorConfig* GetPrecalculatorConfig();

		/*!
		 * @brief Get the Total size of all of the DataSets
		 */
		int GetTotalDataSetSize();

		/*!
		 * @brief Get the size of each of the DataSets
		 */
		vector<int> GetAllDataSetSizes();

		/*!
		 * @brief Get the Start Entries of Each DataSet in the XML
		 */
		vector<int> GetAllStartEntries();

		/*!
		 * @brief Get the Number of Repeats as defined in the top level of the XML
		 */
		int GetNumberRepeats();

		/*!
		 * @brief Get the Random Number Generator Seed as defined in the top level of the XML
		 */
		unsigned int GetSeed();

		/*!
		 * @brief Set the Random Number Generator Seed externally to override any number in the XML
		 */
		void SetSeed( unsigned int new_seed );


		unsigned int GetOriginalSeed() const;
	private:
		/*!
		 * Don't Copy the class this way!
		 */
		XMLConfigReader ( const XMLConfigReader& );

		/*!
		 * Don't Copy the class this way!
		 */
		XMLConfigReader& operator = ( const XMLConfigReader& );

		/*!
		 * All of the file - preciding whitespaces
		 */
		vector<string> wholeFile;

		/*!
		 * All of the tags in the file after it was processed
		 */
		XMLTag* All_XML_Tags;

		/*!
		 * @brief Internal Function to get the ParameterSet from the XML file
		 */
		ParameterSet* GetRawFitParameters();

		PhaseSpaceBoundary* FindCommonPhaseSpace( XMLTag* PhaseTag = NULL );	/*!	Find the Common PhaseSpace for the whole fit, which is top level under RapidFit	*/

		//IPDF* FindCommonPDF( PhaseSpaceBoundary* override );	/*!	Find the CommonPDF and construct it (not currently used due to complex issues!)	*/

		XMLTag* FindCommonPDFXML();		/*!	Finds the CommonPDF object in the XML file which has to be a top level under RapidFit	*/

		XMLTag* FindCommonPDFConfiguratorXML(); /*!	Finds the CommonPDFConfigrator XML object in the file which is top level under RapidFit	*/

		vector< XMLTag* > children;		/*!	All of the Children of the tag wholeFile	*/

		int seed;				/*!	Random Seed	*/

		bool TestXML();				/*!	Internal method for determininig if an XML is valid on class Construction	*/

		string fileName;			/*!	Name of the XML File that this was constructed with	*/
		vector<XMLTag*> fileTags;		/*!	XMLTags extracted from the XML File used as input	*/
		bool XMLValid;				/*! 	Internal boolean to know if the XML is valid or not	*/

		/*!
		 * @brief This is an internal function which appends the label of each PDF to have a unique number which makes deciphering which PDF in a large fit exactly threw the error
		 */
		void AppendCommonNames( XMLTag* input, unsigned int &subParamNum, unsigned int &appendParamNum, unsigned int &configParamNum );

		int _original_seed;
};

#endif


