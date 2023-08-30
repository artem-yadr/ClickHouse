#pragma once

#include <atomic>
#include <variant>
#include <Core/Block.h>
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/multi_polygon.hpp>

#include <Dictionaries/DictionaryStructure.h>
#include <Dictionaries/IDictionary.h>
#include <Dictionaries/IDictionarySource.h>
#include <Dictionaries/DictionaryHelpers.h>


namespace DB
{

namespace bg = boost::geometry;

class IPointDictionary : public IDictionary
{
public:
        enum class InputType
    {
        Array,
        Tuple
    };
    /** Controls the different types allowed for providing the coordinates of points.
      * Right now a point can be represented by either an array or a tuple of two Float64 values.
      */
    enum class PolygonType
    {
        MultiPolygon,
        SimplePolygon
    };

    struct Configuration
    {
        InputType input_type = InputType::Tuple;

        PolygonType polygon_type = PolygonType::SimplePolygon;

        /// Store polygon key column. That will allow to read columns from polygon dictionary.
        bool store_polygon_key_column = false;
    };

    IPointDictionary(
            const StorageID & dict_id_,
            const DictionaryStructure & dict_struct_,
            DictionarySourcePtr source_ptr_,
            DictionaryLifetime dict_lifetime_,
            Configuration configuration_);

    std::string getTypeName() const override { return "Point"; }

    size_t getBytesAllocated() const override { return bytes_allocated; }

    size_t getQueryCount() const override { return query_count.load(std::memory_order_relaxed); }

    double getFoundRate() const override
    {
        size_t queries = query_count.load(std::memory_order_relaxed);
        if (!queries)
            return 0;
        return static_cast<double>(found_count.load(std::memory_order_relaxed)) / queries;
    }

    double getHitRate() const override { return 1.0; }

    size_t getElementCount() const override { return attributes_columns.empty() ? 0 : attributes_columns.front()->size(); }

    double getLoadFactor() const override { return 1.0; }

    DictionarySourcePtr getSource() const override { return source_ptr; }

    const DictionaryStructure & getStructure() const override { return dict_struct; }

    const DictionaryLifetime  & getLifetime() const override { return dict_lifetime; }

    bool isInjective(const std::string & attribute_name) const override { return dict_struct.getAttribute(attribute_name).injective; }

    DictionaryKeyType getKeyType() const override { return DictionaryKeyType::Complex; }

    void convertKeyColumns(Columns & key_columns, DataTypes & key_types) const override;

    ColumnPtr getColumn(
        const std::string& attribute_name,
        const DataTypePtr & result_type,
        const Columns & key_columns,
        const DataTypes & key_types,
        const ColumnPtr & default_values_column) const override;

    ColumnUInt8::Ptr hasKeys(const Columns & key_columns, const DataTypes & key_types) const override;

    Pipe read(const Names & column_names, size_t max_block_size, size_t num_streams) const override;

    /** Single coordinate type. */
    using Coord = Float32;
    /** A two-dimensional point in Euclidean coordinates. */
    using Point = bg::model::d2::point_xy<Coord, bg::cs::cartesian>;
    /** A polygon in boost is a an outer ring of points with zero or more cut out inner rings. */
    using Polygon = bg::model::polygon<Point>;
    /** A ring in boost used for describing the polygons. */
    using Ring = bg::model::ring<Point>;

protected:
    /** Returns true if the given point can be found in the polygon dictionary.
     *  If true id is set to the index of a polygon containing the given point.
     *  Overridden in different implementations of this interface.
     */
    virtual bool find(const Polygon & polygon, size_t & point_index) const = 0;

    std::vector<Point> points;

    const DictionaryStructure dict_struct;
    const DictionarySourcePtr source_ptr;
    const DictionaryLifetime dict_lifetime;
    const Configuration configuration;

private:
/** Helper functions for loading the data from the configuration.
     *  The polygons serving as keys are extracted into boost types.
     *  All other values are stored in one column per attribute.
     */
    void setup();
    void blockToAttributes(const Block & block);
    void loadData();

    void calculateBytesAllocated();

    /** Checks whether a given attribute exists and returns its index */
    size_t getAttributeIndex(const std::string & attribute_name) const;

    /** Helper function for retrieving the value of an attribute by key. */
    template <typename AttributeType, typename ValueGetter, typename ValueSetter, typename DefaultValueExtractor>
    void getItemsImpl(
        const std::vector<IPointDictionary::Point> & requested_key_points,
        ValueGetter && get_value,
        ValueSetter && set_value,
        DefaultValueExtractor & default_value_extractor) const;

    ColumnPtr key_attribute_column;

    Columns attributes_columns;

    size_t bytes_allocated = 0;
    mutable std::atomic<size_t> query_count{0};
    mutable std::atomic<size_t> found_count{0};

    /** Extracts a list of polygons from a column according to input_type and point_type.
     *  The polygons are appended to the dictionary with the corresponding ids.
     */
    void extractPolygons(const ColumnPtr & column);

    /** Extracts a list of points from two columns representing their x and y coordinates. */
    void extractPoints(const ColumnPtr & column);
    
};

}
