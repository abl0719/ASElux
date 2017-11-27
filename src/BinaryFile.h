#ifndef BINARYFILE_H_
#define BINARYFILE_H_

#include <iostream>
#include <fstream>
#include <string>
#include <type_traits>
#include <vector>
#include <cstdint> // if Offset is int32_t instead of ios::streamoff
using namespace std;

template<typename T> struct False : false_type {};

template<typename T>
struct is_string
{
    enum { value = is_same<string, typename remove_cv<T>::type>::value };
};

class BinaryFile
{
public:
    
    typedef ios::streamoff Offset;
    
    ~BinaryFile()
    {
        m_stream.close();
    }
    
    bool openExisting(const std::string& filename)
    {
        if (m_stream.is_open())
            return false;
        m_stream.open(filename, ios::in|ios::out|ios::binary);
        return m_stream.good();
    }
    
    bool createNew(const std::string& filename)
    {
        if (m_stream.is_open())
            return false;
        m_stream.open(filename, ios::in|ios::out|ios::binary|ios::trunc);
        return m_stream.good();
    }
    
    void close()
    {
        if (m_stream.is_open())
            m_stream.close();
    }
    
    template<typename T>
    bool write(const T& data, Offset toOffset)
    {
        static_assert(!is_pointer<T>::value && !is_member_pointer<T>::value,
                      "BinaryFile::write can not be used to write a pointer");
        static_assert(!is_string<T>::value,
                      "BinaryFile::write can not be used to write a std::string");
        static_assert(is_trivially_copyable<T>::value ||
                      is_string<T>::value,  // suppress msg for string
                      "BinaryFile::write can not be used to write a non-trivially copyable class");
        
        return write(reinterpret_cast<const char*>(&data), sizeof(data), toOffset);
    }
    
    template<typename T>
    bool write_vector(vector<T> const& data, Offset toOffset){
        if(data.empty()){
            return false;
        }
        return write(reinterpret_cast<const char*>(data.data()), sizeof(T)*data.size(), toOffset);
    }
    
    template<typename T>
    bool write_array(const T * data, size_t n, Offset toOffset){
        return write(reinterpret_cast<const char*>(data), sizeof(T)*n, toOffset);
    }
    
    bool write(const char* data, size_t length, Offset toOffset)
    {
        return m_stream.seekp(toOffset, ios::beg)  &&
        m_stream.write(data, length);
    }
    
    template<typename T>
    bool read(T& data, Offset fromOffset)
    {
        static_assert(!is_pointer<T>::value && !is_member_pointer<T>::value,
                      "BinaryFile::read can not be used to read a pointer");
        static_assert(!is_string<T>::value,
                      "BinaryFile::read can not be used to read a std::string");
        static_assert(is_trivially_copyable<T>::value ||
                      is_string<T>::value,  // suppress msg for string
                      "BinaryFile::read can not be used to read a non-trivially copyable class");
        
        return read(reinterpret_cast<char*>(&data), sizeof(data), fromOffset);
    }
    
    template<typename T>
    bool read_vector(vector<T> & data, Offset fromOffset){
        return read(reinterpret_cast<char*>(data.data()), sizeof(T)*data.size(), fromOffset);
    }
    
    template<typename T>
    bool read_array(T * data, size_t n, Offset fromOffset){
        return read(reinterpret_cast<char*>(data), sizeof(T)*n, fromOffset);
    }
    
    bool read(char* data, size_t length, Offset fromOffset)
    {
        bool result = m_stream.seekg(fromOffset, ios::beg)  &&
        m_stream.read(data, length);
        if ( ! result)
            m_stream.clear();
        return result;
    }
    
    template<typename T>
    bool read(const T&, Offset)
    {
        static_assert(False<T>::value,
                      "BinaryFile::read can not be used to read into a const or a temporary object");
        return false;
    }
    
    template<typename T>
    bool read(T*&&, Offset)
    {
        static_assert(False<T>::value,
                      "BinaryFile::read can not be used to read a pointer");
        return false;
    }
    
    Offset fileLength()
    {
        if ( ! m_stream.is_open())
            return -1;
        ios::streamoff currPos = m_stream.tellg();
        m_stream.seekg(0, ios::end);
        ios::streamoff length = m_stream.tellg();
        m_stream.seekg(currPos, ios::beg);
        return static_cast<Offset>(length);
    }
    
    bool isOpen() const
    {
        return m_stream.is_open();
    }
    
private:
    fstream m_stream;
    
    // fstreams are not copyable, so BinaryFiles won't be copyable.
};

#endif // BINARYFILE_H_
