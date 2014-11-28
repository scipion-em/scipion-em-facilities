# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'jmdelarosa@cnb.csic.es'
# *
# **************************************************************************
"""
This modules holds the base classes for the ORM implementation.
The Object class is the root in the hierarchy and some other
basic classes.
"""

from itertools import izip

# Binary relations always involve two objects, we 
# call them parent-child objects, the following
# constants reflect which direction of the relation we refer
RELATION_CHILDS = 0
RELATION_PARENTS = 1


class Object(object):
    """ All objects in our Domain should inherit from this class
    that will contains all base properties"""
    def __init__(self, value=None, **args):
        object.__init__(self)
        self._objIsPointer = args.get('objIsPointer', False) # True if will be treated as a reference for storage
        self._objId = args.get('objId', None) # Unique identifier of this object in some context
        self._objParentId = args.get('objParentId', None) # identifier of the parent object
        self._objName = args.get('objName', '') # The name of the object will contains the whole path of ancestors
        self._objLabel = args.get('objLabel', '') # This will serve to label the objects
        self._objComment = args.get('objComment', '')
        self._objTag = args.get('objTag', None) # This attribute serve to make some annotation on the object.
        self._objDoStore = args.get('objDoStore', True) # True if this object will be stored from his parent
        self._objCreation = None
        self._objParent = None # Reference to parent object
        self._objEnabled = True
        self.set(value)

    def getClassName(self):
        return self.__class__.__name__
    
    def getClass(self):
        return type(self)
    
    def getDoc(self):
        return self.__doc__ or ''
    
    def hasAttribute(self, attrName):
        return hasattr(self, attrName)
    
    def getAttributeValue(self, attrName, defaultValue=None):
        """ Get the attribute value given its name.
        Equivalent to getattr(self, name).get() 
        """
        attr = getattr(self, attrName, None)
        if attr is None:
            value = defaultValue
        elif callable(attr):
            value = attr()
        else:
            value = attr.get()
        return value
    
    def setAttributeValue(self, attrName, value):
        """ Set the attribute value given its name.
        Equivalent to setattr(self, name).set(value) 
        If the attrName contains dot: x.y
        it will be equivalent to getattr(getattr(self, 'x'), 'y').set(value)
        """
        attrList = attrName.split('.')
        obj = self
        for attrName in attrList:
            obj = getattr(obj, attrName)
        obj.set(value)
        
    def getAttributes(self):
        """Return the list of attributes than are
        subclasses of Object"""
        for key, attr in self.__dict__.iteritems():
            if issubclass(attr.__class__, Object):
                yield (key, attr)        
                
    def getAttributesToStore(self):
        """Return the list of attributes than are
        subclasses of Object and will be stored"""
        for key, attr in self.getAttributes():
            if attr is not None and attr._objDoStore:
                yield (key, attr)
                
    def isPointer(self):
        """If this is true, the value field is a pointer 
        to another object"""
        return self._objIsPointer
        
    def _convertValue(self, value):
        """Convert a value to desired scalar type"""
        return value
    
    def set(self, value):
        """Set the internal value, if it is different from None
        call the convert function in subclasses"""
        if not value is None:
            value = self._convertValue(value)            
        self._objValue = value
    
    def get(self):
        """Return internal value"""
        return self._objValue
    
    def trace(self, callback):
        """ Add an observer when the set method is called. """
        if self.set == self.__setTrace:
            pass #print "trace already set"
        else:
            self.__set = self.set 
            self.set = self.__setTrace
        self.__setCallback = callback 
        
    def __setTrace(self, value):
        self.__set(value)
        self.__setCallback()
    
    def getObjValue(self):
        """Return the internal value for storage.
        This is a good place to do some update of the
        internal value before been stored"""
        return self._objValue
    
    def getObjId(self):
        """Return object id"""
        return self._objId
    
    def setObjId(self, newId):
        """Set the object id"""
        self._objId = newId
        
    def copyObjId(self, other):
        """ Copy the object id form other to self. """
        self.setObjId(other.getObjId())
        
    def hasObjId(self):
        return not self._objId is None
    
    def cleanObjId(self):
        """ This function will set to None this object id
        and the id of all its children attributes.
        This function should be used when retrieving
        an object from a mapper and storing in a different one.
        """
        self.setObjId(None)
        for _, attr in self.getAttributesToStore():
            attr.cleanObjId()
            
    def getObjParentId(self):
        return self._objParentId
    
    def hasObjParentId(self):
        return self._objParentId is not None
            
    def getObjLabel(self):
        """ Return the label associated with this object"""
        return self._objLabel
    
    def setObjLabel(self, label):
        """ Set the label to better identify this object"""
        self._objLabel = label
             
    def getObjComment(self):
        """ Return the comment associated with this object"""
        return self._objComment
    
    def setObjComment(self, comment):
        """ Set the comment to better identify this object"""
        self._objComment = comment       
        
    def setObjCreation(self, creation):
        """ Set the creation time of the object. """
        self._objCreation = creation
        
    def getObjCreation(self):
        """ Return the stored creation time of the object. """
        return self._objCreation
    
    def strId(self):
        """String representation of id"""
        return str(self._objId)
    
    def getName(self):
        #TODO: REMOVE THIS FUNCTION, SINCE IT DOES NOT COMPLAIN WITH _objX naming
        return self._objName
    
    def getObjName(self):
        return self._objName
    
    def setEnabled(self, enabled):
        self._objEnabled = bool(enabled)
        
    def isEnabled(self):
        """Return if object is enabled"""
        return self._objEnabled
    
    def getNameId(self):
        """ Return an unique and readable id that identifies this object. """
        label = self.getObjLabel()
        if len(label) > 0:
            return label
        elif self.hasObjId():
            return '%s.%s' % (self.getName(), self.strId())
        return ''
    
    def getLastName(self):
        """ If the name contains parent path, remove it
        and only return the attribute name in its parent. 
        """
        if '.' in self._objName:
            return self._objName.split('.')[-1]
        return self._objName 
    
    def setName(self, name):
        self._objName = name
        
    def hasValue(self):        
        return True
    
    def getStore(self):
        """Return True if the object will be stored by the mapper"""
        return self._objDoStore
    
    def setStore(self, value):
        """set the store flag"""
        self._objDoStore = value
    
    def __eq__(self, other):
        """Comparison for scalars should be by value
        and for other objects by reference"""
        if self._objValue is None:
            return object.__eq__(other)
        return self._objValue == other._objValue
    
    def equalAttributes(self, other, ignore=[], verbose=False):
        """Compare that all attributes are equal"""
        for k, v1 in self.getAttributes():
            #v1 = getattr(self, k) # This is necessary because of FakedObject simulation of getattr
            # Skip comparison of attribute names in 'ignore' list
            if k in ignore:
                continue
            v2 = getattr(other, k)
            if issubclass(type(v1), Object):
                comp = v1.equalAttributes(v2, ignore=ignore, verbose=verbose)
            else:
                comp = v1 == v2
            if not comp:
                if verbose:
                    print "Different attributes: "
                    print "self.%s = %s" % (k, v1)
                    print "other.%s = %s" % (k, v2)
                return False
        return True
            
    def copyAttributes(self, other, *attrNames):
        """ Copy attributes in attrNames from other to self. 
        If the name X is in attrNames, it would be equivalent to:
        self.X.set(other.X.get())
        """
        for name in attrNames:
            if isinstance(getattr(other, name), PointerList):
                for pointer in getattr(other, name):
                    getattr(self, name).append(pointer)
            else:
                getattr(self, name).set(getattr(other, name).get())
            
    def __getObjDict(self, prefix, objDict, includeClass):
        if prefix:
            prefix += '.'
        for k, v in self.getAttributesToStore():
            if not v.isPointer():
                kPrefix = prefix + k
                if includeClass:
                    objDict[kPrefix] = (v.getClassName(), v.getObjValue())
                else:
                    objDict[kPrefix] = v.getObjValue()
                if not isinstance(v, Scalar):
                    v.__getObjDict(kPrefix, objDict, includeClass)
            
    def getObjDict(self, includeClass=False):
        """ Return all attributes and values in a dictionary.
        Nested attributes will be separated with a dot in the dict key.
        """
        from collections import OrderedDict
        d = OrderedDict()
        if includeClass:
            d['self'] = (self.getClassName(),)
        self.__getObjDict('', d, includeClass)
        return d
    
    def copy(self, other, copyId=True, ignoreAttrs=[]):
        """ Copy all attributes values from one object to the other.
        The attributes will be created if needed with the corresponding type.
        Params:
            other: the other object from which to make the copy.
            copyId: if true, the _objId will be also copied.
            ignoreAttrs: pass a list with attributes names to ignore.
        """
        copyDict = {'internalPointers': []} 
        self._copy(other, copyDict, copyId, ignoreAttrs=ignoreAttrs)
        self._updatePointers(copyDict)
        return copyDict
        
    def _updatePointers(self, copyDict):
        """ Update the internal pointers after a copy. 
        If there are pointers to other object in the copy 
        the references should be updated.
        """
        for ptr in copyDict['internalPointers']:
            pointedId = ptr.getObjValue().getObjId()
            if  pointedId in copyDict:
                ptr.set(copyDict[pointedId])
        
    def _copy(self, other, copyDict, copyId, level=1, ignoreAttrs=[]):
        """ This method will recursively clone all attributes from one object to the other.
        NOTE: Currently, we are not deleting attributes missing in the 'other' object.
        copyDict: this dict is used to store the ids map between 'other' and 'self' attributes
            This copyDict is used for update pointers and relations later on.
            This will only work if the ids of 'other' attributes has been properly set.
        """
        # Copy basic object data
        #self._objName = other._objName
        if copyId:
            self._objId = other._objId
        self._objValue = other._objValue
        self._objLabel = other._objLabel
        self._objComment = other._objComment
        # Copy attributes recursively
        for name, attr in other.getAttributes():
            if name not in ignoreAttrs:
                myAttr = getattr(self, name, None)
    
                if myAttr is None:
                    myAttr = attr.getClass()()
                    setattr(self, name, myAttr)
                    
                myAttr._copy(attr, copyDict, copyId, level+2)
                # Store the attr in the copyDict
                if attr.hasObjId():
                    #" storing in copyDict with id=", attr.getObjId()
                    copyDict[attr.getObjId()] = myAttr
                # Use the copyDict to fix the reference in the copying object
                # if the pointed one is inside the same object
                if myAttr.isPointer() and myAttr.hasValue():
                    copyDict['internalPointers'].append(myAttr)
    
    def clone(self):
        clone = self.getClass()()
        clone.copy(self)        
        return clone    
    
    def evalCondition(self, condition):
        """ Check if condition is meet.
        Params:
            condition: the condition string, it can contains variables
                or methods without arguments to be evaluated.
            Examples:
                hasCTF
                hasCTF and not hasAligment
        Return:
            The value of the condition evaluated with values
        """
        # Split in possible tokens
        import re
        tokens = re.split('\W+', condition)
        condStr = condition 
        
        for t in tokens:
            if self.hasAttribute(t):
                condStr = condStr.replace(t, str(self.getAttributeValue(t)))
        return eval(condStr)
    
    def printAll(self, name=None, level=0):
        """Print object and all its attributes.
        Mainly for debugging"""
        tab = ' ' * (level*3)
        idStr = ' (id = %s, pid = %s)' % (self.getObjId(), self._objParentId)
        if name is None:
            print tab, self.getClassName(), idStr
        else:
            if name == 'submitTemplate': # Skip this because very large value
                value = '...'
            else:
                value = self.getObjValue()
                
            print tab, '%s = %s' % (name, value), idStr
        for k, v in self.getAttributes():
            v.printAll(k, level + 1)
            
    def printObjDict(self, includeClasses=False):
        """Print object dictionary. Main for debugging"""
        import pprint
        pp = pprint.PrettyPrinter(indent=4)
        pp.pprint(dict(self.getObjDict(includeClasses)))        


class OrderedObject(Object):
    """This is based on Object, but keep the list
    of the attributes to store in the same order
    of insertion, this can be useful where order matters"""
    def __init__(self, value=None, **args):
        object.__setattr__(self, '_attributes', [])
        Object.__init__(self, value, **args)

    def __attrPointed(self, name, value):
        """ Check if a value is already pointed by other
        attribute. This will prevent to storing pointed
        attributes such as:
        self.inputMics = self.inputMicrographs.get()
        In this case we want to avoid to store self.inputMics as 
        another attribute of this object.
        """
        for key in self._attributes:
            attr = getattr(self, key)
            if attr.isPointer():
                if attr.get() is value:
                    return True
        return False
    
    def __setattr__(self, name, value):
        if (not name in self._attributes and
            issubclass(value.__class__, Object) and
            not self.__attrPointed(name, value) and value._objDoStore):
            self._attributes.append(name)
        Object.__setattr__(self, name, value)
    
    def getAttributes(self):
        """Return the list of attributes than are
        subclasses of Object and will be stored"""
        for key in self._attributes:
            yield (key, getattr(self, key))
                
    def deleteAttribute(self, attrName):
        """ Delete an attribute. """
        if attrName in self._attributes:
            self._attributes.remove(attrName)
            delattr(self, attrName)
            
            
class FakedObject(Object):
    """This is based on Object, but will hide the set and get
    access to the attributes, they need to be defined with addAttribute"""
    def __init__(self, value=None, **args):
        object.__setattr__(self, '_attributes', {})
        Object.__init__(self, value, **args)
        
    def addAttribute(self, name, attrClass, **args):
        self._attributes[name] = attrClass(**args)
           
    def __setattr__(self, name, value):
        if name in self._attributes:
            if issubclass(type(value), Object):
                self._attributes[name] = value
            else:
                self._attributes[name].set(value)
        else:
            object.__setattr__(self, name, value)
    
    def __getattr__(self, name):
        if name in self._attributes:
            attr = self._attributes[name]
            if issubclass(type(attr), Scalar):
                return attr.get()
            else:
                return attr
        return None

    def getAttributes(self):
        """Return the list of attributes than are
        subclasses of Object and will be stored"""
        return self._attributes.iteritems()

                
class Scalar(Object):
    """Base class for basic types"""
    def hasValue(self):        
        return self._objValue is not None
    
    def equalAttributes(self, other, ignore=[], verbose=False):
        """Compare that all attributes are equal"""
        return self._objValue == other._objValue
    
    def __str__(self):
        """String representation of the scalar value"""
        return str(self._objValue)
    
    def __eq__(self, other):
        """Comparison for scalars should be by value
        and for other objects by reference"""
        if isinstance(other, Object):
            return self._objValue == other._objValue
        return self._objValue == other

    def __ne__(self, other):
        return not self.__eq__(other)
    
    def __cmp__(self, other):
        """ Comparison implementation for scalars. """
        if isinstance(other, Object):
            return cmp(self._objValue, other._objValue)
        return cmp(self._objValue, other)        
       
    def get(self, default=None):
        """Get the value, if internal value is None
        the default argument passed is returned"""
        if self.hasValue():
            return self._objValue
        return default
    
    def _copy(self, other, *args):
        self.set(other.get())
        
    def swap(self, other):
        """ Swap the contained value between
        self and other objects.
        """
        tmp = self._objValue
        self._objValue = other._objValue
        other._objValue = tmp

    def sum(self, value):
        self._objValue += self._convertValue(value)
        
    def multiply(self, value):
        self._objValue *= value
        
    
class Integer(Scalar):
    """Integer object"""
    def _convertValue(self, value):
        return int(value)
    
    def increment(self):
        """ Add 1 to the current value. """
        self._objValue += 1
    
        
class String(Scalar):
    """String object"""
    def _convertValue(self, value):
        return str(value)
    
    def empty(self):
        """ Return true if None or len == 0 """
        if not self.hasValue():
            return True
        return len(self.get().strip()) == 0
    
        
class Float(Scalar):
    """Float object"""
    EQUAL_PRECISION = 0.001
    
    @classmethod
    def setPrecision(cls, newPrecision):
        """ Set the precision to compare float values.
        Mainly used for testing purposes.
        """
        cls.EQUAL_PRECISION = newPrecision
        
    def _convertValue(self, value):
        return float(value)
    
    def equalAttributes(self, other, ignore=[], verbose=False):
        """Compare that all attributes are equal"""
        # If both float has some value distinct of None
        # then we should compare the absolute value of difference 
        # agains the EQUAL_PRECISION to say equal or not
        if self.hasValue() and other.hasValue():
            return abs(self._objValue - other._objValue) < self.EQUAL_PRECISION
        # If one has None as value, then both should have None
        # to have equal attributes
        if not self.hasValue() and not other.hasValue():
            return True
        
        return False
        
        
class Boolean(Scalar):
    """Boolean object"""
    
    def _convertValue(self, value):
        t = type(value)
        if t is bool:
            return value
        if t is str or t is unicode:
            v = value.strip().lower()
            return v == 'true' or v == '1'
        return bool(value) 
    
    def __nonzero__(self):
        if not self.hasValue():
            return False
        return self.get() 
    
    def __bool__(self):
        return self.get()  
    
    
class Pointer(Object):
    """Reference object to other one"""
    def __init__(self, value=None, **args):
        Object.__init__(self, value, objIsPointer=True, **args)
        # The _extended attribute will be used to point to attributes of a pointed object
        # or the id of an item inside a set
        self._extended = String() 
       
    def __str__(self):
        """String representation of a pointer"""
        if self.hasValue():
            className = self.getObjValue().getClassName()
            strId = self.getObjValue().strId()
            return '-> %s (%s)' % (className, strId)
        return '-> None'

    def hasValue(self):
        return self._objValue is not None
    
    def get(self, default=None):
        """ Get the pointed object. 
        If the _extended property is set, then the extended attribute or item id
        will be retrieved by the call to get().
        """
        extended = self._extended.get()
        if extended:
            if extended.startswith('__attribute__'):
                attribute = extended.split('__')[-1]
                value = getattr(self._objValue, attribute, default)
            elif extended.startswith('__itemid__'):
                itemId = int(extended.split('__')[-1])
                value = self._objValue[itemId]
                value._parentObject = self._objValue
            else:
                # This is now only for compatibility reasons
                try:
                    itemId = int(extended)
                    value = self._objValue[itemId]
                    value._parentObject = self._objValue
                except Exception:
                    attribute = str(extended)
                    value = getattr(self._objValue, attribute, default)
                    
                #raise Exception("Invalid value '%s' for pointer._extended property." % extended)
        else:
            value = self._objValue
            
        return value
    
    def set(self, other):
        """ Set the pointer value but cleanning the extendend property. """
        Object.set(self, other)
        # This check is needed because set is call from the Object constructor
        # when this attribute is not setup yet (a dirty patch, I know)
        if hasattr(self, '_extended'):
            self._extended.set(None)
        
    def setExtendedAttribute(self, attributeName):
        """ Point to an attribute of the pointed object. """
        self._extended.set('__attribute__' + attributeName)
        
    def setExtendedItemId(self, itemId):
        """ Point to an specific item of a pointed Set. """
        self._extended.set('__itemid__%d' % itemId)
        
    def getAttributes(self):
        yield ('_extended', getattr(self, '_extended'))
    

class List(Object, list):
    ITEM_PREFIX = '__item__'
    
    """Class to store a list of objects"""
    def __init__(self, **args):
        Object.__init__(self, **args)
        list.__init__(self)
        
    def __getattr__(self, name):
        if name.startswith(self.ITEM_PREFIX):
            i = self._stringToIndex(name)
            if i < len(self):
                return self[i]
        raise AttributeError("List object has not attribute: " + name)
            
    def __setattr__(self, name, value):
        if name.startswith('__item__') or len(name)==0:
            self.append(value)
        else:
            object.__setattr__(self, name, value)

    def getAttributes(self):
        # First yield all attributes not contained in the list
        for name, attr in Object.getAttributes(self):
            yield (name, attr)
        # Now yield elements contained in the list
        for i, item in enumerate(self):
            yield (self._indexToString(i), item)
            
    def _indexToString(self, i):
        """Return the way the string index is generated.
        String indexes will start in 1, that's why i+1
        """
        return "%s%06d" % (self.ITEM_PREFIX, i+1)
    
    def _stringToIndex(self, strIndex):
        """ From the string index representation obtain the index.
        For simetry the number in the index string will be 
        decreased in 1.
        """
        return int(strIndex.split(self.ITEM_PREFIX)[1]) - 1
            
    #TODO: check if needed
    def __len__(self):
        return list.__len__(self)
    
    def isEmpty(self):
        return len(self) > 0
    
    def clear(self):
        del self[:]
        
    def _convertValue(self, value):
        """Value should be a list."""
        if not isinstance(value, list):
            raise Exception("List.set: value should be a list.")
        self.clear()
        for item in value:
            self.append(item)
        return None
    
        
class PointerList(List):
    def __init__(self, **kwargs):
        List.__init__(self, **kwargs)
        
    def _convertValue(self, value):
        if isinstance(value, list):
            self.clear()
            for obj in value:
                self.append(Pointer(value=obj))
        else:
            raise Exception("Could not set a PointerList value to: %s" % value)

    def append(self, value):
        """ Append Pointer of objects to the list.
         If value is a Pointer, just add it to the list.
         If is another subclass of Object, create
         a Pointer first and append the pointer.
        """
        if isinstance(value, Pointer):
            pointer = value
        elif isinstance(value, Object):
            pointer = Pointer()
            pointer.set(value)
        else:
            raise Exception("Only subclasses of Object can be added to PointerList")

        List.append(self, pointer)

            
class CsvList(Scalar, list):
    """This class will store a list of objects
    in a single DB row separated by comma.
    pType: the type of the list elememnts, int, bool, str"""
    def __init__(self, pType=str, **kwargs):
        Scalar.__init__(self, **kwargs)
        list.__init__(self)
        self._pType = pType
        
    def _convertValue(self, value):
        """Value should be a str with comman separated values
        or a list.
        """
        self.clear()
        if value:
            if isinstance(value, str) or isinstance(value, unicode):
                for s in value.split(','):
                    self.append(self._pType(s))
            elif isinstance(value, list) or isinstance(value, tuple):
                for s in value:
                    self.append(self._pType(s))
            else:
                raise Exception("CsvList.set: Invalid value type: ", type(value))
            
    def getObjValue(self):
        self._objValue = ','.join(map(str, self))
        return self._objValue
    
    def get(self):
        return self.getObjValue()
    
    def __str__(self):
        return list.__str__(self)
     
    def isEmpty(self):
        return len(self) == 0
    
    def clear(self):
        del self[:]
        
        
class Array(Object):
    """Class for holding fixed len array"""
    def __init__(self, size=10, **args):
        Object.__init__(self, size, **args)
        
    def set(self, size):
        """Set the array size"""
        self._objValue = int(size)  
        for i in range(int(size)):
            self.__setitem__(i, None)                 
        
    def strIndex(self, i):
        return 'item_%04d' % i
    
    def __setitem__(self, index, value):
        self.__dict__[self.strIndex(index)] = value
        
    def __getitem__(self, index):
        return self.__dict__[self.strIndex(index)]
    
    def __len__(self):
        return self._objValue
    
    
class Set(OrderedObject):
    """ This class will be a container implementation for elements.
    It will use an extra sqlite file to store the elements.
    All items will have an unique id that identifies each element in the set.
    """
    ITEM_TYPE = None # This property should be defined to know the item type
    
    def __init__(self, filename=None, prefix='', 
                 mapperClass=None, classesDict=None, **args):
        # Use the object value to store the filename
        OrderedObject.__init__(self, **args)
        self._mapper = None
        self._idCount = 0
        self._size = Integer(0) # cached value of the number of images  
        #self._idMap = {}#FIXME, remove this after id is the one in mapper
        self.setMapperClass(mapperClass)
        self._mapperPath = CsvList() # sqlite filename
        self._mapperPath.trace(self.load) # Load the mapper whenever the filename is changed
        self._representative = None
        self._classesDict = classesDict 
        # If filename is passed in the constructor, it means that
        # we want to create a new object, so we need to delete it if
        # the file exists
        if filename:
            self._mapperPath.set('%s, %s' % (filename, prefix)) # This will cause the creation of the mapper

    def aggregate(self, operations
                      , operationLabel
                      , groupByLabels=None):
        return self._mapper.aggregate(operations, operationLabel, groupByLabels)

    def setMapperClass(self, MapperClass):
        """ Set the mapper to be used for storage. """
        if MapperClass is None:
            from pyworkflow.mapper.sqlite import SqliteFlatMapper
            MapperClass = SqliteFlatMapper
        Object.__setattr__(self, '_MapperClass', MapperClass)
        
    def __getitem__(self, itemId):
        """ Get the image with the given id. """
        return self._mapper.selectById(itemId)

    def __contains__(self, itemId):
        """ element in Set """
        return self._mapper.selectById(itemId) != None

    def _iterItems(self, random=False):
        return self._mapper.selectAll(random=random)#has flat mapper, iterate is true

    def getFirstItem(self):
        """ Return the first item in the Set. """
        return self._mapper.selectFirst()
    
    def __iter__(self):
        """ Iterate over the set of images. """
        return self._iterItems()
       
    def __len__(self):
        return self._size.get()
    
    def getSize(self):
        """Return the number of images"""
        return self._size.get()
    
    def getFileName(self):
        if len(self._mapperPath):
            return self._mapperPath[0]
        return None
    
    def getPrefix(self):
        if len(self._mapperPath) > 1:
            return self._mapperPath[1]
        return None
    
    def write(self):
        """This method will be used to persist in a file the
        list of images path contained in this Set
        path: output file path
        images: list with the images path to be stored
        """
        #TODO: If mapper is in memory, do commit and dump to disk
        self._mapper.setProperty('self', self.getClassName())
        objDict = self.getObjDict()
        for key, value in objDict.iteritems():
            self._mapper.setProperty(key, value)
        self._mapper.commit()
    
    def _loadClassesDict(self):
        return self._classesDict or globals()
    
    def setClassesDict(self, classesDict):
        """ Set the dictionary with classes where to look for classes names. """
        self._classesDict = classesDict
    
    def load(self):
        """ Load extra data from files. """
        if self._mapperPath.isEmpty():
            raise Exception("Set.load:  mapper path and prefix not set.")
        fn, prefix = self._mapperPath
        self._mapper = self._MapperClass(fn, self._loadClassesDict(), prefix)
        self._size.set(self._mapper.count())
           
    def close(self):
        self._mapper.close()
        
    def clear(self):
        self._mapper.clear()
        self._idCount = 0
        self._size.set(0)
         
    def append(self, item):
        """ Add an item to the set.
        If the item has already an id, use it.
        If not, keep a counter with the max id
        and assign the next one.
        """
        # The _idCount and _size properties work fine
        # under the assumption that once a Set is stored,
        # then it is read-only (no more appends).
        #
        # Anyway, this can be easily changed by updating
        # both from the underlying sqlite when reading a set.

        if not item.hasObjId():
            self._idCount += 1
            item.setObjId(self._idCount)
        else:
            self._idCount = max(self._idCount, item.getObjId()) + 1
        self._insertItem(item)
        self._size.increment()
#        self._idMap[item.getObjId()] = item
        
    def _insertItem(self, item):
        self._mapper.insert(item)
        
    def update(self, item):
        """ Update an existing item. """
        self._mapper.update(item)
                
    def __str__(self):
        return "%-20s (%d items)" % (self.getClassName(), self.getSize())
    
    def getSubset(self, n):
        """ Return a subset of n element, making a clone of each. """
        subset = []
        for i, item in enumerate(self):
            subset.append(item.clone())
            if i == n:
                break
        return subset
    
    def setRepresentative(self, representative):
        self._representative = representative
    
    def getRepresentative(self):       
        return self._representative
    
    def hasRepresentative(self):
        """ Return true if have a representative image. """
        return self._representative is not None

    def equalItemAttributes(self, other, ignore=[], verbose=False):
        """Compare that all items in self and other
        return True for equalAttributes.
        """
        return all(x.getObjId() == y.getObjId() and
                   x.equalAttributes(y, ignore=ignore, verbose=verbose)
                   for x, y in izip(self, other))


def ObjectWrap(value):
    """This function will act as a simple Factory
    to create objects from Python basic types"""
    t = type(value)
    if issubclass(t, Object):
        return value
    if t is int or t is long:
        return Integer(value)
    if t is bool:
        return Boolean(value)
    if t is float:
        return Float(value)
    if t is list:
        o = CsvList()
        o.set(value)
        return o
    if t is None:
        return None
    #If it is str, unicode or unknown type, convert to string
    return String(value)
         
           