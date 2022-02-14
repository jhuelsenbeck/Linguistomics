#include "container.hpp"


BufferAllocator::BufferAllocator(size_t size, size_t elements) {
    startBuffer = nullptr;
    endBuffer = nullptr;
    elementSize = size;
    numElements = 0;
    bufSize = 0;
    create(elements);
}

BufferAllocator::BufferAllocator(const BufferAllocator &b) {
    startBuffer = nullptr;
    endBuffer = nullptr;
    elementSize = b.elementSize;
    numElements = 0;
    bufSize = 0;
    copy(b);
}

BufferAllocator::~BufferAllocator() {
    deallocate();
}

void BufferAllocator::deallocate() {
    if (startBuffer) 
        {
        free(startBuffer);
        startBuffer = nullptr;
        endBuffer   = nullptr;
        numElements = 0;
        bufSize = 0;
        }
}

void BufferAllocator::create(size_t elements) {
    if (elements != numElements)
    {
        if (numElements > 0)
            deallocate();
        if (elements > 0) {
            bufSize = elementSize * elements;
            startBuffer = (char*)malloc(bufSize);
        }
        numElements = elements;
        endBuffer = startBuffer + bufSize;
    }
}

void BufferAllocator::setZero() {
    memset(startBuffer, 0, bufSize);
}

void BufferAllocator::copy(const BufferAllocator &b) {

    create(b.numElements);
    memcpy(startBuffer, b.startBuffer, bufSize);
}

bool BufferAllocator::operator==(const BufferAllocator& b) const {
    if (numElements != b.numElements)
        return false;
    return memcmp(startBuffer, b.startBuffer, bufSize) == 0;
}

bool BufferAllocator::operator!=(const BufferAllocator& b) const {

    if (numElements != b.numElements)
        return true;
    return memcmp(startBuffer, b.startBuffer, bufSize) != 0;
}


