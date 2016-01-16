/****************************************************************************
** Meta object code from reading C++ file 'renderthread.h'
**
** Created: Sat Jan 16 12:36:20 2016
**      by: The Qt Meta Object Compiler version 63 (Qt 4.8.4)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "renderthread.h"
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'renderthread.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 63
#error "This file was generated using the moc from 4.8.4. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
static const uint qt_meta_data_RenderThread[] = {

 // content:
       6,       // revision
       0,       // classname
       0,    0, // classinfo
       2,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       2,       // signalCount

 // signals: signature, parameters, type, tag, flags
      32,   14,   13,   13, 0x05,
      61,   13,   13,   13, 0x05,

       0        // eod
};

static const char qt_meta_stringdata_RenderThread[] = {
    "RenderThread\0\0image,scaleFactor\0"
    "renderedImage(QImage,double)\0stopped()\0"
};

void RenderThread::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        Q_ASSERT(staticMetaObject.cast(_o));
        RenderThread *_t = static_cast<RenderThread *>(_o);
        switch (_id) {
        case 0: _t->renderedImage((*reinterpret_cast< const QImage(*)>(_a[1])),(*reinterpret_cast< double(*)>(_a[2]))); break;
        case 1: _t->stopped(); break;
        default: ;
        }
    }
}

const QMetaObjectExtraData RenderThread::staticMetaObjectExtraData = {
    0,  qt_static_metacall 
};

const QMetaObject RenderThread::staticMetaObject = {
    { &QThread::staticMetaObject, qt_meta_stringdata_RenderThread,
      qt_meta_data_RenderThread, &staticMetaObjectExtraData }
};

#ifdef Q_NO_DATA_RELOCATION
const QMetaObject &RenderThread::getStaticMetaObject() { return staticMetaObject; }
#endif //Q_NO_DATA_RELOCATION

const QMetaObject *RenderThread::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->metaObject : &staticMetaObject;
}

void *RenderThread::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_RenderThread))
        return static_cast<void*>(const_cast< RenderThread*>(this));
    return QThread::qt_metacast(_clname);
}

int RenderThread::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QThread::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 2)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 2;
    }
    return _id;
}

// SIGNAL 0
void RenderThread::renderedImage(const QImage & _t1, double _t2)
{
    void *_a[] = { 0, const_cast<void*>(reinterpret_cast<const void*>(&_t1)), const_cast<void*>(reinterpret_cast<const void*>(&_t2)) };
    QMetaObject::activate(this, &staticMetaObject, 0, _a);
}

// SIGNAL 1
void RenderThread::stopped()
{
    QMetaObject::activate(this, &staticMetaObject, 1, 0);
}
QT_END_MOC_NAMESPACE
