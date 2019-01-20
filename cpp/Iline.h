// header file for declarations of Ifile class

#ifndef ILINE
#define ILINE

#include <iostream>
#include <string>
#include <vector>

//=====================================================================================================
/**
 * Class to hold information and manipulate text lines. Reunite tools to parse files by make tests with
 * the text lines content and features.
 * @class Iline
 * @author Igor Barden Grillo
 * @date 07/04/18
 * @file Iline.h
 * @brief class to represent and manipulate a line from text.
 */
class Iline {
	public:
		//---------------------------------------------------------------
		/**
		 * @brief STL vector to hold the individual words of the line.
		 */
		std::vector<std::string> words;
		
		//---------------------------------------------------------------
		/**
		 * @brief Integer with the line length.
		 */
		unsigned int line_len;
		
		//---------------------------------------------------------------
		/**
		 * Instantiates an empty object adn initializes the member data to default values.
		 * @brief Default constructor of Iline class.
		 */
		Iline();
		
		//---------------------------------------------------------------
		/**
		 * Intantiates the object from a string with the line iformation to be stored and manipulated.
		 * The line is tokenized and the words is stored in words member data STL vector of strings.
		 * @brief Constructor from a string with the line content.
		 * @param String reference with line content.
		 */		
		Iline(std::string Line);
		
		//---------------------------------------------------------------
		/**
		 * Instatiates the object from copying the member data of another class type object.
		 * @brief Copy constructor of the class.
		 * @param Ilne Const reference from the object to be copied.
		 */
		Iline(const Iline& rhs_line);
		
		//---------------------------------------------------------------
		/**
		 * Return this object with member data overwritted with member data from another 
		 * class type object passed in the right hand side of assigment operator.
		 * @brief Assigment operator overloading 
		 * @param Iline Const reference from the object to be copied.
		 * @return Ilne const reference from this object.
		 */
		Iline& operator=(const Iline& rhs_line); 
		
		/**
		 * Instantiates an Iline object from stolen ownership of memeber data 
		 * from a rvalue reference of an class type object.
		 * @brief Move constructor of the class. 
		 * @param Iline rvalue reference from the object to be moved.
		 */
		Iline(Iline&& rhs_line) noexcept;
		 
		/**
		 * Return this object with member data overwrittd from stolen memeber data ownership
		 * of another class type object.
		 * @brief Move assigment operator overloading.
		 * @param Iline rvalue reference from the object to be moved.
		 * @return Ilne const reference from this object.
		 */
		Iline& operator=(Iline&& rhs_line) noexcept;
		
		//---------------------------------------------------------------
		/**
		 * Memeber function to test if an word, or a sequence of caracteres with whitespaces before
		 * and until the position passed in the function is present on the line
		 * @brief Class method to test if the string and the position passed is the word in the line. 
		 * @param String to be tested.
		 * @param Integer with the position of the word in the line.
		 * @param Integer value to sinalize the end of the word to be tested in the line. 
		 * @return Boolean variable.
		 */		 
		bool IF_word(std::string s,int pos,int fin);
		
		//---------------------------------------------------------------
		/**
		 * Overloaded member function to test if the entire line matchs the string passe to the function. 
		 * @brief Test if the line have the same content of the string passed.
		 * @param String to test if its equal to the line in the current object. 
		 * @return Boolean variable indicating the result of the test.
		 */ 
		bool IF_line(std::string s);
		
		//---------------------------------------------------------------
		/**
		 * Member function overloaded to test if the line has one of its tokenized words
		 * in the position passed in the arguments.
		 * @brief Test if the line have the word in the position passed. 
		 * @param String to test if is in the line.
		 * @param Integer value with the indice of the position of the word in the line.
		 * @param Lenght of the line to compare with the cuurrent object.
		 * @return Boolean variable indicating the result of the test.
		 */
		bool IF_line(std::string s, int pos,int len);
		
		//---------------------------------------------------------------
		/**
		 * Member function overloaded to test if the line has two of its tokenized words
		 * in the positions passed in the arguments.
		 * @brief Test if the line have the words in the positions passed. 
		 * @param String to test if is in the line.
		 * @param Integer value with the indice of the first word in the line.
		 * @param String to test if is in the line.
		 * @param Integer value with the indice of the second word in the line.
		 * @param Lenght of the line to compare with the cuurrent object.
		 * @return Boolean variable.
		 */
		bool IF_line(std::string s1, int pos1,std::string s2, int pos2, int len);
		
		//---------------------------------------------------------------
		/**
		 * @brief Return the line by combining the elements of the vector words.
		 * @return String with the words of the line. 
		 */
		std::string get_line(); 
		
		//---------------------------------------------------------------
		/**
		 * This member function erases the element coverted of the words vector and realocates 
		 * all elements to new positions.
		 * @brief Convert a string in words vector to integer, returns and erases from the conteiner.
		 * @param Integer with the position of the element in the vector.
		 * @return Integer value from the converted vector element.
		 */
		int pop_int(int pos);
		
		//---------------------------------------------------------------
		/**
		 * This member function does not erases the element in the vector.
		 * @brief Convert a string in words vector to integer and return from the conteiner.
		 * @param Integer with the position of the element in the vector.
		 * @return Integer value from the converted vector element.
		 */
		int get_int(int pos);
		
		//---------------------------------------------------------------
		/**
		 * This member function erases the element coverted of the words vector and realocates 
		 * all elements to new positions.
		 * @brief Convert a string in words vector to double, return and erases from the conteiner.
		 * @param Integer with the position of the element in the vector.
		 * @return Double value from the converted vector element.
		 */
		double pop_double(int pos);
		
		//---------------------------------------------------------------
		/**
		 * This member function does not erases the element in the vector. 
		 * @brief Convert a string in words vector to double, return and erases from the conteiner.
		 * @param Integer with the position of the element in the vector.
		 * @return Double with the converted value.
		 */
		double get_double(int pos);
		
		//---------------------------------------------------------------
		/**
		 * This member function erases the element coverted of the words vector and realocates 
		 * all elements to new positions.
		 * @brief Convert a from words vector string to double with double precision signalized with a "D" 
		 * instead a "E" and erase the element from the words vector.
		 * @param Integer with the position of the element in the vector.
		 * @return Double with the converted value.
		 */
		double pop_double_f(int pos);
		
		//---------------------------------------------------------------
		/**
		 * This member function does not erases the element in the vector.
		 * @brief Convert a from words vector string to double with double precision signalized with a "D" 
		 * instead a "E".
		 * @param Integer with the position of the element in the vector.
		 * @return Double with the converted value.
		 */
		double get_double_f(int pos);
		
		//---------------------------------------------------------------
		/**
		 * This member function erases the element coverted of the words vector and realocates 
		 * all elements to new positions.
		 * @brief 
		 * @param Integer with the position of the element in the vector.
		 * @return 
		 */
		std::string pop_string(int pos);
		
		//---------------------------------------------------------------
		/**
		 * This member function does not erases the element in the vector.
		 * @brief Returns a c
		 * @param Integer with the position of the element in the vector.
		 * @return String from the words vector.
		 */
		std::string& get_string(int pos);
				
		//---------------------------------------------------------------
		/**
		 * @brief Class method to print to the console line information.
		 * @return None.
		 */
		void print();
		
		//---------------------------------------------------------------
		/**
		 * Destructor of the class.
		 */	
		~Iline();
	
};

#endif