#pragma once
#include <ostream>

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
