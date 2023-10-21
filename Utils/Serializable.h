#ifndef SERIALIZABLE_H
#define SERIALIZABLE_H

#include "Epsilon.h"
#include <Eigen/Dense>
#include <iostream>
#include <memory>

namespace PHYSICSMOTION {
struct IOData;
struct SerializableBase {
  explicit SerializableBase() {}
  virtual ~SerializableBase() {}
  virtual bool read(std::istream& is);
  virtual bool write(std::ostream& os) const;
  virtual bool read(std::istream& is,IOData* dat);
  virtual bool write(std::ostream& os,IOData* dat) const;
  virtual bool readStr(const std::string& path);
  virtual bool writeStr(const std::string& path) const;
  virtual std::shared_ptr<SerializableBase> copy() const;
  virtual std::string type() const=0;
  virtual bool operator<(const SerializableBase& other) const {
    return type()<other.type();
  }
};

//minimal serializable system
std::shared_ptr<IOData> getIOData();
std::ostream& writeSerializableData(std::shared_ptr<SerializableBase> val,std::ostream& os,IOData* dat=NULL);
std::istream& readSerializableData(std::shared_ptr<SerializableBase>& val,std::istream& is,IOData* dat=NULL);
void registerType(IOData* dat,std::shared_ptr<SerializableBase> type);
void registerTypeAs(IOData* dat,std::shared_ptr<SerializableBase> alias,std::shared_ptr<SerializableBase> type);
template <typename T>inline void registerType(IOData* dat) {
  registerType(dat,std::shared_ptr<SerializableBase>(new T));
}
template <typename A,typename T>inline void registerTypeAs(IOData* dat) {
  registerTypeAs(dat,std::shared_ptr<SerializableBase>(new A),std::shared_ptr<SerializableBase>(new T));
}

//io for shared_ptr
template <typename T>inline std::ostream& writeBinaryData(std::shared_ptr<T> val,std::ostream& os,IOData* dat=NULL) {
  std::shared_ptr<SerializableBase> valS=std::dynamic_pointer_cast<SerializableBase>(val);
  return writeSerializableData(valS,os,dat);
}
template <typename T>inline std::istream& readBinaryData(std::shared_ptr<T>& val,std::istream& is,IOData* dat=NULL) {
  std::shared_ptr<SerializableBase> valS;
  readSerializableData(valS,is,dat);
  val=std::dynamic_pointer_cast<T>(valS);
  return is;
}
template <typename T>inline std::ostream& writeBinaryData(std::weak_ptr<T> val,std::ostream& os,IOData* dat=NULL) {
  std::shared_ptr<T> val_shared=val.lock();
  std::shared_ptr<SerializableBase> valS=std::dynamic_pointer_cast<SerializableBase>(val_shared);
  return writeSerializableData(valS,os,dat);
}
template <typename T>inline std::istream& readBinaryData(std::weak_ptr<T>& val,std::istream& is,IOData* dat=NULL) {
  std::shared_ptr<SerializableBase> valS;
  readSerializableData(valS,is,dat);
  val=std::dynamic_pointer_cast<T>(valS);
  return is;
}
extern std::ostream& writeBinaryData(const SerializableBase& val,std::ostream& os,IOData* dat=NULL);
extern std::istream& readBinaryData(SerializableBase& val,std::istream& is,IOData* dat=NULL);
}

#endif
