/**
 * @file Manager-inl.h
 * @author  Scott Rasmussen (scott.rasmussen@zaita.com)
 * @version 1.0
 * @date 22/10/2012
 * @section LICENSE
 *
 * Copyright NIWA Science �2012 - www.niwa.co.nz
 *
 * $Date: 2008-03-04 16:33:32 +1300 (Tue, 04 Mar 2008) $
 */

// Namespaces
namespace isam {
namespace base {

/**
 * This method is our singleton instance method
 */
template <class ClassType, class StoredType>
ClassType& Manager<ClassType, StoredType>::Instance() {
  static ClassType singleton;
  return singleton;
}

/**
 * Validate the objects stored in the manager
 */
template <class ClassType, class StoredType>
void Manager<ClassType, StoredType>::Validate() {
  for(shared_ptr<StoredType> stored_object : objects_) {
    stored_object->Validate();
  }
}

/**
 * Build the objects stored in the manager
 */
template <class ClassType, class StoredType>
void Manager<ClassType, StoredType>::Build() {
  for(shared_ptr<StoredType> stored_object : objects_) {
    stored_object->Build();
  }
}

} /* namespae base */
} /* namespace isam */
