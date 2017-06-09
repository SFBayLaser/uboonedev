///////////////////////////////////////////////////////////////////////
///
/// \file   IDraw2DObjects.h
///
/// \brief  This provides an interface for tools which are tasked with
///         drawing objects in the 2D event display
///
/// \author T. Usher
///
////////////////////////////////////////////////////////////////////////

#ifndef IDraw2DObjects_H
#define IDraw2DObjects_H

#include "art/Framework/Principal/Event.h"

namespace display
{
    class IDraw2DObjects
    {
    public:
        virtual ~IDraw2DObjects() noexcept = default;
        
        // Basic drawing of recob::Wire object
        virtual void drawWire2D(const art::Event&) const = 0;
    };
}

#endif
