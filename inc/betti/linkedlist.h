#ifndef LINKED_LIST_H
#define LINKED_LIST_H

/**
 * Node class for double linked list
 */
template <class T>
class LL_Node {
public:
	LL_Node(T _obj) : next(NULL), prev(NULL) {
		object = _obj;
	}

	LL_Node<T>* prev;
	LL_Node<T>* next;
	T	object;
};

/**
 * LinkedList -- Double Linked List
 */
template <class T>
class LinkedList {
public:
	/**
	 * Construct an empty list
	 */
	LinkedList();

	/**
	 * Destructor
	 */
	~LinkedList();

	/**
	 * Return the head of the list
	 * @return NULL if the list is empty
	 */
	LL_Node<T>* head() const {
		return p_head;
	}

	/**
	 * Return the tail of the list
	 * @return NULL if the list is empty
	 */
	LL_Node<T>* tail() const {
		return p_tail;
	}

	/**
	 * Insert a node pointer at the end of the list.
	 * @note Linked list will take the responsibility of releasing the node
	 */
	void insert(LL_Node<T>* p_node);

	/**
	 * Remove the node pointed by the given pointer.
	 * @note Memory will be released.
	 */
	void remove(LL_Node<T>* p_node);

	/**
	 * Check if the list is empty.
	 */
	bool empty() const {
		return (p_head == NULL);
	}

	/**
	 * Check if a node pointer is the tail
	 */
	bool isTail(LL_Node<T>* p_node) {
		return (p_node == p_tail);
	}

	/**
	 */
	int size() const {
		return m_size;
	}

private:
	LL_Node<T>* p_head;
	LL_Node<T>* p_tail;
	int m_size;
};

/*
 *	Double linked list
 */
template <class T>
LinkedList<T>::LinkedList()
{
	p_head = NULL;
	p_tail = NULL;
	m_size = 0;
}

template <class T>
LinkedList<T>::~LinkedList()
{
	if(!empty()) {
		LL_Node<T>* p_node = head();
		while(p_node != NULL) {
			LL_Node<T>* tmp_ptr = p_node->next;
			delete p_node;
			p_node = tmp_ptr;
		}
	}
}

template <class T>
void LinkedList<T>::insert(LL_Node<T>* p_node)
{
	p_node->next = NULL;
	p_node->prev = p_tail;
	if(p_tail != NULL) p_tail->next = p_node;
	if(p_head == NULL) p_head = p_node;
	p_tail = p_node;
	m_size ++;
}

template <class T>
void LinkedList<T>::remove(LL_Node<T>* p_node)
{
	if(p_node->prev == NULL) { // head
		p_head = p_node->next;
	} else {
		p_node->prev->next = p_node->next;
	}

	if(p_node->next == NULL) { // tail
		p_tail = p_node->prev;
	} else {
		p_node->next->prev = p_node->prev;
	}
	m_size --;
	delete p_node;
}
#endif



