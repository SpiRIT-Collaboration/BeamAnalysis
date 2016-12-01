// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME dOlibbeam_dict

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "RConfig.h"
#include "TClass.h"
#include "TDictAttributeMap.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TInterpreter.h"
#include "TVirtualMutex.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TIsAProxy.h"
#include "TFileMergeInfo.h"
#include <algorithm>
#include "TCollectionProxyInfo.h"
/*******************************************************************/

#include "TDataMember.h"

// Since CINT ignores the std namespace, we need to do so in this file.
namespace std {} using namespace std;

// Header files passed as explicit arguments
#include "BeamRaw.h"
#include "BeamBeam.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static TClass *BeamRaw_Dictionary();
   static void BeamRaw_TClassManip(TClass*);
   static void *new_BeamRaw(void *p = 0);
   static void *newArray_BeamRaw(Long_t size, void *p);
   static void delete_BeamRaw(void *p);
   static void deleteArray_BeamRaw(void *p);
   static void destruct_BeamRaw(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::BeamRaw*)
   {
      ::BeamRaw *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::BeamRaw));
      static ::ROOT::TGenericClassInfo 
         instance("BeamRaw", "BeamRaw.h", 13,
                  typeid(::BeamRaw), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &BeamRaw_Dictionary, isa_proxy, 0,
                  sizeof(::BeamRaw) );
      instance.SetNew(&new_BeamRaw);
      instance.SetNewArray(&newArray_BeamRaw);
      instance.SetDelete(&delete_BeamRaw);
      instance.SetDeleteArray(&deleteArray_BeamRaw);
      instance.SetDestructor(&destruct_BeamRaw);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::BeamRaw*)
   {
      return GenerateInitInstanceLocal((::BeamRaw*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::BeamRaw*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *BeamRaw_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::BeamRaw*)0x0)->GetClass();
      BeamRaw_TClassManip(theClass);
   return theClass;
   }

   static void BeamRaw_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *BeamBeam_Dictionary();
   static void BeamBeam_TClassManip(TClass*);
   static void *new_BeamBeam(void *p = 0);
   static void *newArray_BeamBeam(Long_t size, void *p);
   static void delete_BeamBeam(void *p);
   static void deleteArray_BeamBeam(void *p);
   static void destruct_BeamBeam(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::BeamBeam*)
   {
      ::BeamBeam *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::BeamBeam));
      static ::ROOT::TGenericClassInfo 
         instance("BeamBeam", "BeamBeam.h", 13,
                  typeid(::BeamBeam), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &BeamBeam_Dictionary, isa_proxy, 0,
                  sizeof(::BeamBeam) );
      instance.SetNew(&new_BeamBeam);
      instance.SetNewArray(&newArray_BeamBeam);
      instance.SetDelete(&delete_BeamBeam);
      instance.SetDeleteArray(&deleteArray_BeamBeam);
      instance.SetDestructor(&destruct_BeamBeam);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::BeamBeam*)
   {
      return GenerateInitInstanceLocal((::BeamBeam*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::BeamBeam*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *BeamBeam_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::BeamBeam*)0x0)->GetClass();
      BeamBeam_TClassManip(theClass);
   return theClass;
   }

   static void BeamBeam_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_BeamRaw(void *p) {
      return  p ? new(p) ::BeamRaw : new ::BeamRaw;
   }
   static void *newArray_BeamRaw(Long_t nElements, void *p) {
      return p ? new(p) ::BeamRaw[nElements] : new ::BeamRaw[nElements];
   }
   // Wrapper around operator delete
   static void delete_BeamRaw(void *p) {
      delete ((::BeamRaw*)p);
   }
   static void deleteArray_BeamRaw(void *p) {
      delete [] ((::BeamRaw*)p);
   }
   static void destruct_BeamRaw(void *p) {
      typedef ::BeamRaw current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::BeamRaw

namespace ROOT {
   // Wrappers around operator new
   static void *new_BeamBeam(void *p) {
      return  p ? new(p) ::BeamBeam : new ::BeamBeam;
   }
   static void *newArray_BeamBeam(Long_t nElements, void *p) {
      return p ? new(p) ::BeamBeam[nElements] : new ::BeamBeam[nElements];
   }
   // Wrapper around operator delete
   static void delete_BeamBeam(void *p) {
      delete ((::BeamBeam*)p);
   }
   static void deleteArray_BeamBeam(void *p) {
      delete [] ((::BeamBeam*)p);
   }
   static void destruct_BeamBeam(void *p) {
      typedef ::BeamBeam current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::BeamBeam

namespace {
  void TriggerDictionaryInitialization_dOlibbeam_dict_Impl() {
    static const char* headers[] = {
"BeamRaw.h",
"BeamBeam.h",
0
    };
    static const char* includePaths[] = {
"/mnt/research/spirit/SPIRIT_TPC/20160919.v1.03/include/root",
"/mnt/research/spirit/SPIRIT_USERS/barneyj/anaroot/BeamAnalysis/BeamAnalysis/macros/",
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "dOlibbeam_dict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
class __attribute__((annotate("$clingAutoload$BeamRaw.h")))  BeamRaw;
class __attribute__((annotate("$clingAutoload$BeamBeam.h")))  BeamBeam;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "dOlibbeam_dict dictionary payload"

#ifndef G__VECTOR_HAS_CLASS_ITERATOR
  #define G__VECTOR_HAS_CLASS_ITERATOR 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
#include "BeamRaw.h"
#include "BeamBeam.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
"BeamBeam", payloadCode, "@",
"BeamRaw", payloadCode, "@",
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule(".libbeam_dict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_dOlibbeam_dict_Impl, {}, classesHeaders);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_dOlibbeam_dict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_dOlibbeam_dict() {
  TriggerDictionaryInitialization_dOlibbeam_dict_Impl();
}
