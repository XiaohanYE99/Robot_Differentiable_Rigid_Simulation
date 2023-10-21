#include "Serializable.h"
#include "Pragma.h"
#include "IO.h"
#include <fstream>
#include <unordered_map>

namespace PHYSICSMOTION {
//SerializableBase
bool SerializableBase::read(std::istream& is) {
  FUNCTION_NOT_IMPLEMENTED
  return false;
}
bool SerializableBase::write(std::ostream& os) const {
  FUNCTION_NOT_IMPLEMENTED
  return false;
}
bool SerializableBase::read(std::istream& is,IOData* dat) {
  return read(is);
}
bool SerializableBase::write(std::ostream& os,IOData* dat) const {
  return write(os);
}
bool SerializableBase::readStr(const std::string& path) {
  std::ifstream is(path,std::ios::binary);
  std::shared_ptr<IOData> dat=getIOData();
  return read(is,dat.get());
}
bool SerializableBase::writeStr(const std::string& path) const {
  std::ofstream os(path,std::ios::binary);
  std::shared_ptr<IOData> dat=getIOData();
  return write(os,dat.get());
}
std::shared_ptr<SerializableBase> SerializableBase::copy() const {
  FUNCTION_NOT_IMPLEMENTED
  return std::shared_ptr<SerializableBase>();
}

//IOData
struct IOData {
 public:
  typedef std::shared_ptr<SerializableBase> SPTR;
  typedef std::allocator<std::pair<const SPTR,int> > ALLOC_MAP_WR;
  typedef std::allocator<std::pair<const int,SPTR> > ALLOC_MAP_RD;
  typedef std::allocator<std::pair<const std::string,SPTR> > ALLOC_TYPE_SET;
  typedef std::unordered_map<SPTR,int,std::hash<SPTR>,std::equal_to<SPTR>,ALLOC_MAP_WR> MAP_WR;
  typedef std::unordered_map<int,SPTR,std::hash<int>,std::equal_to<int>,ALLOC_MAP_RD> MAP_RD;
  typedef std::unordered_map<std::string,SPTR,std::hash<std::string>,std::equal_to<std::string>,ALLOC_TYPE_SET> TYPE_SET;
 public:
  IOData():_index(0) {}
  int getIndex() {
    return _index++;
  }
  void registerType(std::shared_ptr<SerializableBase> type) {
    ASSERT_MSG(!type->type().empty(),"Given type doesn't support shared_ptr serialization!");
    if(_tSet.find(type->type()) != _tSet.end()) {
      ASSERT_MSGV(typeid(*_tSet[type->type()]) == typeid(*type),"Conflicit type id: %s",type->type().c_str())
    }
    _tSet[type->type()]=type;
  }
  void registerTypeAs(std::shared_ptr<SerializableBase> alias,std::shared_ptr<SerializableBase> type) {
    ASSERT_MSG(!alias->type().empty(),"Given type doesn't support shared_ptr serialization!");
    if(_tSet.find(alias->type()) != _tSet.end()) {
      ASSERT_MSGV(typeid(*_tSet[alias->type()]) == typeid(*type),"Conflicit type id: %s",alias->type().c_str())
    }
    _tSet[alias->type()]=type;
  }
  template <typename T>void createNew(std::istream& is,std::shared_ptr<T>& val) const {
    std::string type;
    readBinaryData(type,is);
    ASSERT_MSG(type != "","Type not found!")
    for(TYPE_SET::const_iterator beg=_tSet.begin(),end=_tSet.end(); beg!=end; beg++)
      if(beg->first == type) {
        val=std::dynamic_pointer_cast<T>(beg->second->copy());
        return;
      }
    ASSERT_MSG(false,"Cannot find compatible type!")
  }
  MAP_WR _ptrMapWr;
  MAP_RD _ptrMapRd;
 private:
  int _index;
  TYPE_SET _tSet;
};
//minimal serializable system
std::shared_ptr<IOData> getIOData() {
  return std::shared_ptr<IOData>(new IOData);
}
std::ostream& writeSerializableData(std::shared_ptr<SerializableBase> val,std::ostream& os,IOData* dat) {
  ASSERT_MSG(dat,"You must provide pool for serialize shared_ptr!")
  if(val.get() == NULL) {
    writeBinaryData((int)-1,os);
    return os;
  }
  std::shared_ptr<SerializableBase> ptrS=std::dynamic_pointer_cast<SerializableBase>(val);
  ASSERT_MSG(ptrS,"Not serializable type!")
  IOData::MAP_WR::const_iterator it=dat->_ptrMapWr.find(ptrS);
  if(it == dat->_ptrMapWr.end()) {
    int id=dat->getIndex();
    writeBinaryData(id,os,dat);
    writeBinaryData(val->type(),os,dat);
    dat->_ptrMapWr[val]=id;
    writeBinaryData(*val,os,dat);
  } else {
    writeBinaryData(it->second,os,dat);
  }
  return os;
}
std::istream& readSerializableData(std::shared_ptr<SerializableBase>& val,std::istream& is,IOData* dat) {
  ASSERT_MSG(dat,"You must provide pool for serialize shared_ptr!")
  int id;
  readBinaryData(id,is,dat);
  if(id == (int)-1) {
    val.reset((SerializableBase*)NULL);
    return is;
  }
  IOData::MAP_RD::const_iterator it=dat->_ptrMapRd.find(id);
  if(it == dat->_ptrMapRd.end()) {
    dat->createNew(is,val);
    dat->_ptrMapRd[id]=val;
    readBinaryData(*val,is,dat);
  } else {
    val=std::dynamic_pointer_cast<SerializableBase>(dat->_ptrMapRd[id]);
  }
  return is;
}
void registerType(IOData* dat,std::shared_ptr<SerializableBase> type) {
  dat->registerType(type);
}
void registerTypeAs(IOData* dat,std::shared_ptr<SerializableBase> alias,std::shared_ptr<SerializableBase> type) {
  dat->registerTypeAs(alias,type);
}
std::ostream& writeBinaryData(const SerializableBase& val,std::ostream& os,IOData* dat) {
  val.write(os,dat);
  return os;
}
std::istream& readBinaryData(SerializableBase& val,std::istream& is,IOData* dat) {
  val.read(is,dat);
  return is;
}
}
