#pragma once
#include <ostream>

/**
 * A very simple JSON writer
 * 
 * Assumes UTF-8 encoded strings.
 * Add values with value(), null(), or using the << operator.
 * For arrays, begin with bArray() and end with eArray().
 * For objects, begin with bObject(), set key with safeKey() (known values, no escape) or unsafeKey() (escape) and values with value(), and end with eObject().
 * 
 * The output is written in the stream directly and no state is kept.
 * This class doesn't keep track of nesting and will never error on incorrect use.
 */
class JsonWriter {
    std::ostream &os;
    bool first = true;

    JsonWriter& item();
public:
    explicit JsonWriter(std::ostream &os);

    JsonWriter& value(float v);
    JsonWriter& value(bool v);
    JsonWriter& value(uint32_t v);
    JsonWriter& value(int32_t v);
    JsonWriter& value(const char* v);
    JsonWriter& null();

    JsonWriter& safeKey(const char* v);
    JsonWriter& unsafeKey(const char* v);
    JsonWriter& bObject();
    JsonWriter& eObject();
    JsonWriter& bArray();
    JsonWriter& eArray();
};

template<typename T>
JsonWriter& operator<<(JsonWriter& j, const T& v) {
    return j.value(v);
}
