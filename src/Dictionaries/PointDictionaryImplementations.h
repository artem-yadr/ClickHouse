#pragma once

#include "PointDictionary.h"

#include <vector>

namespace DB
{

/** Simple implementation of the polygon dictionary. Doesn't generate anything during its construction.
  * Iterates over all stored polygons for each query, checking each of them in linear time.
  * Retrieves the polygon with the smallest area containing the given point.
  * If there is more than one any such polygon may be returned.
  */
class PointDictionarySimple : public IPointDictionary
{
public:
    PointDictionarySimple(
            const StorageID & dict_id_,
            const DictionaryStructure & dict_struct_,
            DictionarySourcePtr source_ptr_,
            DictionaryLifetime dict_lifetime_,
            Configuration configuration_);

    std::shared_ptr<const IExternalLoadable> clone() const override;

private:
    bool find(const Polygon & polygon, size_t & point_index) const override;
};

}

